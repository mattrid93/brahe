#include "brahe/event_detector.h"

#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace brahe {

namespace {

constexpr size_t kMaxChildrenPerBody = 64;

bool is_finite_state(const State2& s) {
    return std::isfinite(s.r.x) && std::isfinite(s.r.y) &&
           std::isfinite(s.v.x) && std::isfinite(s.v.y);
}

// Initialise out_event so all fields are finite and deterministic, regardless of the
// code path taken before the first successful assignment. Critical for the
// "no NaNs on ordinary invalid inputs" contract.
void init_output(PredictedEvent& out, double fallback_time) {
    out.type = EventType::None;
    out.time = fallback_time;
    out.from_body = InvalidBody;
    out.to_body = InvalidBody;
    out.state = State2{{0.0, 0.0}, {0.0, 0.0}};
}

// Propagate the initial state forward by dt. dt == 0 returns the initial state
// unchanged (avoids unnecessary propagator calls at start-time boundary checks).
SolveStatus propagate_by(double mu, const State2& initial, double dt, State2& out) {
    if (dt == 0.0) {
        out = initial;
        return SolveStatus::Ok;
    }
    return TwoBody::propagate(mu, initial, dt, out);
}

}  // namespace

EventDetector::EventDetector(const BodySystem& bodies, const Tolerances& tolerances)
    : bodies_(bodies), tolerances_(tolerances) {}

SolveStatus EventDetector::find_next_event(const EventSearchRequest& req,
                                           PredictedEvent& out) const {
    // --- Step 0: always leave output with finite values ---
    const double fallback_time =
        std::isfinite(req.start_time) ? req.start_time : 0.0;
    init_output(out, fallback_time);

    // --- Step 1: validate ---
    if (!std::isfinite(req.start_time) || !std::isfinite(req.time_limit)) {
        return SolveStatus::InvalidInput;
    }
    if (req.time_limit < req.start_time) {
        return SolveStatus::InvalidInput;
    }
    if (!is_finite_state(req.initial_state)) {
        return SolveStatus::InvalidInput;
    }

    const BodyDef* cb = bodies_.get_body(req.central_body);
    if (cb == nullptr) {
        return SolveStatus::InvalidInput;
    }

    // Zero-radius initial state is invalid for two-body propagation.
    if (length_squared(req.initial_state.r) == 0.0) {
        return SolveStatus::InvalidInput;
    }

    const double mu = cb->mu;
    const double body_radius = cb->radius;
    const double soi_radius = cb->soi_radius;
    const BodyId parent_id = cb->parent_id;
    const double root_eps = tolerances_.root_epsilon;
    const double time_eps = tolerances_.time_epsilon;
    const int max_iter = tolerances_.max_event_refine_iterations;
    const double horizon = req.time_limit - req.start_time;

    // Children of the central body (in canonical BodyId order).
    const BodyId* child_begin = bodies_.children_begin(req.central_body);
    const BodyId* child_end = bodies_.children_end(req.central_body);
    const size_t num_children = static_cast<size_t>(child_end - child_begin);
    if (num_children > kMaxChildrenPerBody) {
        return SolveStatus::InvalidInput;
    }

    // --- Step 2: set default output to TimeLimit at horizon end with propagated state ---
    // Propagator may fail for degenerate inputs (e.g. purely radial trajectories).
    // In that case we still return a valid TimeLimit event using the initial state
    // as the fallback -- the detector contract forbids NaN output on valid requests.
    bool propagation_ok = true;
    {
        State2 at_end;
        SolveStatus ps = propagate_by(mu, req.initial_state, horizon, at_end);
        out.type = EventType::TimeLimit;
        out.time = req.time_limit;
        out.from_body = req.central_body;
        out.to_body = InvalidBody;
        if (ps != SolveStatus::Ok || !is_finite_state(at_end)) {
            out.state = req.initial_state;
            propagation_ok = false;
        } else {
            out.state = at_end;
        }
    }

    // --- Step 3: candidate arbitration ---
    PredictedEvent best{};
    bool have_best = false;
    auto consider = [&](const PredictedEvent& ev) {
        if (!have_best || detail::event_precedes(ev, best, time_eps)) {
            best = ev;
            have_best = true;
        }
    };

    // Root evaluation helpers. These operate on a pre-propagated spacecraft state.
    auto eval_impact = [&](const State2& s) { return length(s.r) - body_radius; };
    auto eval_exit = [&](const State2& s) { return length(s.r) - soi_radius; };
    auto eval_child = [&](const State2& s, BodyId child_id, double abs_t) {
        const BodyDef* ch = bodies_.get_body(child_id);
        Vec2 rc = bodies_.position_in_parent(child_id, abs_t);
        return length(s.r - rc) - ch->soi_radius;
    };

    // --- Step 4: start-state / boundary policy (spec 6.3) ---

    // Impact at start: fire unconditionally if already at-or-inside the body radius.
    {
        double f = eval_impact(req.initial_state);
        if (f <= root_eps) {
            PredictedEvent ev;
            ev.type = EventType::Impact;
            ev.time = req.start_time;
            ev.from_body = req.central_body;
            ev.to_body = InvalidBody;
            ev.state = req.initial_state;
            consider(ev);
        }
    }

    // SOI exit at start: fire only if outside by > root_eps OR on boundary with outward trend.
    {
        double f = eval_exit(req.initial_state);
        double r_dot_v = dot(req.initial_state.r, req.initial_state.v);
        bool outside = f > root_eps;
        bool on_boundary_outward = std::abs(f) <= root_eps && r_dot_v > 0.0;
        if (outside || on_boundary_outward) {
            PredictedEvent ev;
            ev.type = EventType::SoiExit;
            ev.time = req.start_time;
            ev.from_body = req.central_body;
            ev.to_body = parent_id;  // InvalidBody if central is root
            ev.state = req.initial_state;
            consider(ev);
        }
    }

    // Child SOI entry at start: fire only if inside by > root_eps OR on boundary with inward trend.
    for (const BodyId* p = child_begin; p != child_end; ++p) {
        const BodyDef* ch = bodies_.get_body(*p);
        if (ch == nullptr) continue;
        Vec2 r_child = bodies_.position_in_parent(*p, req.start_time);
        Vec2 v_child = bodies_.velocity_in_parent(*p, req.start_time);
        Vec2 d = req.initial_state.r - r_child;
        Vec2 v_rel = req.initial_state.v - v_child;
        double dist = length(d);
        double f = dist - ch->soi_radius;
        double d_dot_v = dot(d, v_rel);  // sign of d/dt(dist^2)/2 = sign of d/dt(dist)
        bool inside = f < -root_eps;
        bool on_boundary_inward = std::abs(f) <= root_eps && d_dot_v < 0.0;
        if (inside || on_boundary_inward) {
            PredictedEvent ev;
            ev.type = EventType::SoiEntry;
            ev.time = req.start_time;
            ev.from_body = req.central_body;
            ev.to_body = *p;
            ev.state = req.initial_state;
            consider(ev);
        }
    }

    // Degenerate horizon: skip scan entirely.
    if (horizon <= time_eps) {
        if (have_best) out = best;
        return SolveStatus::Ok;
    }

    // --- Step 5: coarse scan step size ---
    // Adaptive step bounded by smallest child SOI / max relative speed, plus safety caps.
    double v_sc = length(req.initial_state.v);
    double v_rel_max = v_sc;
    double min_child_soi = std::numeric_limits<double>::infinity();
    for (const BodyId* p = child_begin; p != child_end; ++p) {
        const BodyDef* ch = bodies_.get_body(*p);
        if (ch == nullptr) continue;
        double v_moon = std::abs(ch->orbit_radius * ch->angular_rate);
        v_rel_max = std::max(v_rel_max, v_sc + v_moon);
        min_child_soi = std::min(min_child_soi, ch->soi_radius);
    }
    const double v_floor = 0.1;
    double dt_child = std::isfinite(min_child_soi)
                          ? (min_child_soi / std::max(v_rel_max, v_floor))
                          : std::numeric_limits<double>::infinity();
    double dt_default = horizon / 50.0;
    double dt_floor = std::max(time_eps * 100.0, 1e-12);
    double dt_scan = std::min(dt_child, dt_default);
    dt_scan = std::max(dt_scan, dt_floor);
    if (dt_scan > horizon) dt_scan = horizon;

    // --- Step 6: coarse scan loop ---
    double t_lo = req.start_time;
    State2 s_lo = req.initial_state;
    double f_imp_lo = eval_impact(s_lo);
    double f_exit_lo = eval_exit(s_lo);
    double f_child_lo[kMaxChildrenPerBody];
    {
        size_t i = 0;
        for (const BodyId* p = child_begin; p != child_end; ++p, ++i) {
            f_child_lo[i] = eval_child(s_lo, *p, t_lo);
        }
    }

    // Hard cap on scan iterations prevents infinite loops from pathological inputs.
    const int max_scan_steps = 200000;
    int scan_steps = 0;

    // Note: horizon propagation may have failed (e.g. radial trajectory that
    // reaches the central singularity). The scan can still succeed for earlier
    // t values, and the mid-scan failure handler below will then kick in.
    (void)propagation_ok;

    while (t_lo < req.time_limit && scan_steps < max_scan_steps) {
        ++scan_steps;
        double t_hi = std::min(t_lo + dt_scan, req.time_limit);
        State2 s_hi;
        SolveStatus ps = propagate_by(mu, req.initial_state, t_hi - req.start_time, s_hi);
        if (ps != SolveStatus::Ok || !is_finite_state(s_hi)) {
            // Propagator failed mid-scan. The only failure mode after validation is
            // the radial singularity at r=0 -- the spacecraft falls into the central
            // body during (t_lo, t_hi]. Binary-search for the last propagable time
            // and, if an impact bracket exists there, refine it.
            double t_valid_lo = t_lo;   // propagation works here
            double t_valid_hi = t_hi;   // propagation fails here
            State2 s_last = s_lo;
            for (int k = 0; k < 64; ++k) {
                if (t_valid_hi - t_valid_lo <= time_eps) break;
                double t_m = 0.5 * (t_valid_lo + t_valid_hi);
                State2 s_m;
                SolveStatus ms =
                    propagate_by(mu, req.initial_state, t_m - req.start_time, s_m);
                if (ms == SolveStatus::Ok && is_finite_state(s_m)) {
                    t_valid_lo = t_m;
                    s_last = s_m;
                } else {
                    t_valid_hi = t_m;
                }
            }
            double f_imp_end = eval_impact(s_last);
            if (f_imp_lo > 0.0 && f_imp_end <= 0.0) {
                auto f_of_t = [&](double t) {
                    State2 s;
                    (void)propagate_by(mu, req.initial_state, t - req.start_time, s);
                    return eval_impact(s);
                };
                double t_root = t_lo;
                SolveStatus rs = detail::refine_root_bisection(
                    f_of_t, t_lo, t_valid_lo, f_imp_lo, f_imp_end, root_eps, max_iter,
                    t_root);
                if (rs == SolveStatus::Ok) {
                    State2 s_root;
                    (void)propagate_by(mu, req.initial_state, t_root - req.start_time,
                                       s_root);
                    PredictedEvent ev;
                    ev.type = EventType::Impact;
                    ev.time = t_root;
                    ev.from_body = req.central_body;
                    ev.to_body = InvalidBody;
                    ev.state = s_root;
                    consider(ev);
                }
            }
            if (have_best) out = best;
            return SolveStatus::Ok;
        }

        // --- Impact bracket ---
        double f_imp_hi = eval_impact(s_hi);
        if (f_imp_lo > 0.0 && f_imp_hi <= 0.0) {
            auto f_of_t = [&](double t) {
                State2 s;
                (void)propagate_by(mu, req.initial_state, t - req.start_time, s);
                return eval_impact(s);
            };
            double t_root = t_lo;
            SolveStatus rs = detail::refine_root_bisection(
                f_of_t, t_lo, t_hi, f_imp_lo, f_imp_hi, root_eps, max_iter, t_root);
            if (rs == SolveStatus::Ok) {
                State2 s_root;
                (void)propagate_by(mu, req.initial_state, t_root - req.start_time, s_root);
                PredictedEvent ev;
                ev.type = EventType::Impact;
                ev.time = t_root;
                ev.from_body = req.central_body;
                ev.to_body = InvalidBody;
                ev.state = s_root;
                consider(ev);
            }
            // Refinement failure: leave the default TimeLimit / best-so-far output alone.
        }

        // --- SOI exit bracket ---
        double f_exit_hi = eval_exit(s_hi);
        if (f_exit_lo < 0.0 && f_exit_hi >= 0.0) {
            auto f_of_t = [&](double t) {
                State2 s;
                (void)propagate_by(mu, req.initial_state, t - req.start_time, s);
                return eval_exit(s);
            };
            double t_root = t_lo;
            SolveStatus rs = detail::refine_root_bisection(
                f_of_t, t_lo, t_hi, f_exit_lo, f_exit_hi, root_eps, max_iter, t_root);
            if (rs == SolveStatus::Ok) {
                State2 s_root;
                (void)propagate_by(mu, req.initial_state, t_root - req.start_time, s_root);
                PredictedEvent ev;
                ev.type = EventType::SoiExit;
                ev.time = t_root;
                ev.from_body = req.central_body;
                ev.to_body = parent_id;
                ev.state = s_root;
                consider(ev);
            }
            // Refinement failure: leave the default TimeLimit / best-so-far output alone.
        }

        // --- Child entry brackets ---
        size_t ci = 0;
        for (const BodyId* p = child_begin; p != child_end; ++p, ++ci) {
            BodyId child_id = *p;
            double f_ch_lo = f_child_lo[ci];
            double f_ch_hi = eval_child(s_hi, child_id, t_hi);
            if (f_ch_lo > 0.0 && f_ch_hi <= 0.0) {
                auto f_of_t = [&, child_id](double t) {
                    State2 s;
                    (void)propagate_by(mu, req.initial_state, t - req.start_time, s);
                    return eval_child(s, child_id, t);
                };
                double t_root = t_lo;
                SolveStatus rs = detail::refine_root_bisection(
                    f_of_t, t_lo, t_hi, f_ch_lo, f_ch_hi, root_eps, max_iter, t_root);
                if (rs == SolveStatus::Ok) {
                    State2 s_root;
                    (void)propagate_by(mu, req.initial_state, t_root - req.start_time, s_root);
                    PredictedEvent ev;
                    ev.type = EventType::SoiEntry;
                    ev.time = t_root;
                    ev.from_body = req.central_body;
                    ev.to_body = child_id;
                    ev.state = s_root;
                    consider(ev);
                }
                // Refinement failure: leave the default TimeLimit / best-so-far alone.
            }
            f_child_lo[ci] = f_ch_hi;
        }

        // --- Advance scan window ---
        t_lo = t_hi;
        s_lo = s_hi;
        f_imp_lo = f_imp_hi;
        f_exit_lo = f_exit_hi;

        // Early termination: if we already have a confirmed event earlier than the
        // current scan position, no later bracket can beat it (scan is time-monotonic).
        if (have_best && best.time < t_lo - time_eps) {
            break;
        }
    }

    if (have_best) {
        out = best;
    }
    // Otherwise out stays as TimeLimit at horizon end (set in Step 2).
    return SolveStatus::Ok;
}

// --- detail helpers ---

namespace detail {

bool event_precedes(const PredictedEvent& a, const PredictedEvent& b, double time_epsilon) {
    // If times differ by more than epsilon, earlier wins.
    if (a.time < b.time - time_epsilon) return true;
    if (b.time < a.time - time_epsilon) return false;

    // Within epsilon: type priority (lower enum value = higher priority).
    // Impact(1) > SoiEntry(2) > SoiExit(3) > Burn(4) > TimeLimit(5).
    if (a.type != b.type) return static_cast<int>(a.type) < static_cast<int>(b.type);

    // Same type within epsilon: lowest BodyId wins (for sibling SOI entries).
    return a.to_body < b.to_body;
}

}  // namespace detail

}  // namespace brahe
