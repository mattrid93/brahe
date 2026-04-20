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

struct CachedChild {
    BodyId id = InvalidBody;
    double soi_radius = 0.0;
    double orbit_radius = 0.0;
    double angular_rate = 0.0;
    double phase_at_epoch = 0.0;
};

struct CandidateSet {
    bool have_best = false;
    PredictedEvent best{};

    void consider(const PredictedEvent& ev, double time_epsilon) {
        if (!have_best || detail::event_precedes(ev, best, time_epsilon)) {
            best = ev;
            have_best = true;
        }
    }
};

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

PredictedEvent make_event(EventType type, double time, BodyId from_body,
                          BodyId to_body, const State2& state) {
    PredictedEvent ev;
    ev.type = type;
    ev.time = time;
    ev.from_body = from_body;
    ev.to_body = to_body;
    ev.state = state;
    return ev;
}

class CachedConicPropagator {
public:
    CachedConicPropagator(double mu, const State2& initial)
        : mu_(mu), initial_(initial) {
        const double r_mag = length(initial.r);
        const double v_mag = length(initial.v);
        const double h = TwoBody::specific_angular_momentum_z(initial);
        const double h_threshold = 1e-12 * std::max(1.0, r_mag * v_mag);
        if (std::abs(h) <= h_threshold) {
            return;
        }

        SolveStatus status = TwoBody::to_elements(mu, initial, elements_);
        use_cached_ = status == SolveStatus::Ok &&
                      elements_.semi_latus_rectum > 0.0 &&
                      std::isfinite(elements_.semi_latus_rectum) &&
                      std::isfinite(elements_.eccentricity) &&
                      std::isfinite(elements_.angular_momentum_z);
    }

    SolveStatus propagate(double dt, State2& out) const {
        if (dt == 0.0) {
            out = initial_;
            return SolveStatus::Ok;
        }
        if (!use_cached_) {
            return TwoBody::propagate(mu_, initial_, dt, out);
        }

        ConicElements2D el = elements_;
        const double e = el.eccentricity;
        const double p = el.semi_latus_rectum;
        const double sign_h = (el.angular_momentum_z >= 0.0) ? 1.0 : -1.0;
        double nu_new = 0.0;

        if (el.type == ConicType::Ellipse) {
            const double a = el.semi_major_axis;
            if (a <= 0.0 || !std::isfinite(a)) return SolveStatus::NumericalFailure;
            const double n = std::sqrt(mu_ / (a * a * a));
            const double M_new = el.mean_anomaly + sign_h * n * dt;

            double E_new = 0.0;
            SolveStatus ks = detail::mean_to_eccentric_anomaly(M_new, e, E_new);
            if (ks != SolveStatus::Ok) return ks;
            nu_new = detail::eccentric_to_true_anomaly(E_new, e);
        } else if (el.type == ConicType::Hyperbola) {
            const double a = el.semi_major_axis;
            if (a >= 0.0 || !std::isfinite(a)) return SolveStatus::NumericalFailure;
            const double n = std::sqrt(mu_ / (-a * a * a));
            const double M_new = el.mean_anomaly + sign_h * n * dt;

            double H_new = 0.0;
            SolveStatus ks = detail::mean_to_hyperbolic_anomaly(M_new, e, H_new);
            if (ks != SolveStatus::Ok) return ks;
            nu_new = detail::hyperbolic_to_true_anomaly(H_new, e);
        } else {
            if (p <= 0.0 || !std::isfinite(p)) return SolveStatus::NumericalFailure;
            const double D0 = std::tan(el.true_anomaly / 2.0);
            const double coeff = std::sqrt(2.0 * mu_ / (p * p * p));
            const double W = D0 + D0 * D0 * D0 / 3.0 + sign_h * coeff * dt;
            const double D1 = 2.0 * std::sinh(std::asinh(1.5 * W) / 3.0);
            nu_new = 2.0 * std::atan(D1);
        }

        el.true_anomaly = nu_new;
        return TwoBody::from_elements(el, out);
    }

private:
    double mu_ = 0.0;
    State2 initial_{};
    ConicElements2D elements_{};
    bool use_cached_ = false;
};

PredictedEvent make_time_limit_event(const EventSearchRequest& req,
                                     const CachedConicPropagator& propagator,
                                     double horizon) {
    State2 at_end;
    SolveStatus ps = propagator.propagate(horizon, at_end);
    if (ps != SolveStatus::Ok || !is_finite_state(at_end)) {
        at_end = req.initial_state;
    }
    return make_event(EventType::TimeLimit, req.time_limit, req.central_body,
                      InvalidBody, at_end);
}

Vec2 child_position(const CachedChild& child, double t) {
    double theta = child.phase_at_epoch + child.angular_rate * t;
    return {child.orbit_radius * std::cos(theta), child.orbit_radius * std::sin(theta)};
}

Vec2 child_velocity(const CachedChild& child, double t) {
    double theta = child.phase_at_epoch + child.angular_rate * t;
    double v_mag = child.orbit_radius * child.angular_rate;
    return {v_mag * (-std::sin(theta)), v_mag * std::cos(theta)};
}

struct RadiusCrossing {
    bool found = false;
    double dt = 0.0;
};

struct RadialRange {
    double min_radius = 0.0;
    double max_radius = std::numeric_limits<double>::infinity();
};

struct ScanConfig {
    bool scan_impact = false;
    bool scan_exit = false;
    double spacecraft_speed = 0.0;
    double max_relative_speed = 0.0;
    double dt = 0.0;
};

struct ChildCache {
    CachedChild children[kMaxChildrenPerBody];
    size_t count = 0;
    bool radially_culled = false;
    bool have_radial_range = false;
    RadialRange radial_range;
};

enum class RadiusCrossingDirection {
    Inward,
    Outward,
};

void consider_radius_candidate(const ConicElements2D& el, double target_true_anomaly,
                               RadiusCrossingDirection direction, double horizon,
                               double time_eps, RadiusCrossing& best) {
    const double e = el.eccentricity;
    const double sign_h = (el.angular_momentum_z >= 0.0) ? 1.0 : -1.0;
    const double radial_trend = sign_h * std::sin(target_true_anomaly);
    if (direction == RadiusCrossingDirection::Inward && radial_trend >= 0.0) return;
    if (direction == RadiusCrossingDirection::Outward && radial_trend <= 0.0) return;

    double dt = std::numeric_limits<double>::infinity();

    if (el.type == ConicType::Ellipse) {
        const double a = el.semi_major_axis;
        if (a <= 0.0 || !std::isfinite(a)) return;

        const double E = detail::true_to_eccentric_anomaly(target_true_anomaly, e);
        const double M = detail::eccentric_to_mean_anomaly(E, e);
        const double n = std::sqrt(el.mu / (a * a * a));
        double dM = sign_h * (M - el.mean_anomaly);
        const double min_dM = time_eps * n;
        while (dM <= min_dM) dM += 2.0 * M_PI;
        dt = dM / n;
    } else if (el.type == ConicType::Hyperbola) {
        const double a = el.semi_major_axis;
        if (a >= 0.0 || !std::isfinite(a)) return;

        const double H = detail::true_to_hyperbolic_anomaly(target_true_anomaly, e);
        const double M = detail::hyperbolic_to_mean_anomaly(H, e);
        const double n = std::sqrt(el.mu / (-a * a * a));
        dt = (M - el.mean_anomaly) / (sign_h * n);
    } else {
        const double p = el.semi_latus_rectum;
        if (p <= 0.0 || !std::isfinite(p)) return;

        const double D = std::tan(target_true_anomaly / 2.0);
        const double M = D + D * D * D / 3.0;
        const double coeff = std::sqrt(2.0 * el.mu / (p * p * p));
        dt = (M - el.mean_anomaly) / (sign_h * coeff);
    }

    if (!std::isfinite(dt) || dt <= time_eps || dt > horizon + time_eps) return;
    if (!best.found || dt < best.dt) {
        best.found = true;
        best.dt = std::min(dt, horizon);
    }
}

RadiusCrossing first_central_radius_crossing(double mu, const State2& initial,
                                             double radius,
                                             RadiusCrossingDirection direction,
                                             double horizon, double time_eps) {
    RadiusCrossing best;
    if (mu <= 0.0 || radius <= 0.0 || horizon <= time_eps ||
        !std::isfinite(mu) || !std::isfinite(radius)) {
        return best;
    }

    ConicElements2D el;
    if (TwoBody::to_elements(mu, initial, el) != SolveStatus::Ok) return best;
    if (el.semi_latus_rectum <= 0.0 || !std::isfinite(el.semi_latus_rectum)) return best;

    const double e = el.eccentricity;
    if (e <= 1e-12 || !std::isfinite(e)) return best;

    const double cos_nu = (el.semi_latus_rectum / radius - 1.0) / e;
    if (cos_nu < -1.0 - 1e-12 || cos_nu > 1.0 + 1e-12) return best;

    const double nu_abs = detail::safe_acos(cos_nu);
    consider_radius_candidate(el, nu_abs, direction, horizon, time_eps, best);
    if (nu_abs > 0.0 && nu_abs < M_PI) {
        consider_radius_candidate(el, -nu_abs, direction, horizon, time_eps, best);
    }
    return best;
}

bool conic_radial_range(double mu, const State2& initial, RadialRange& out) {
    ConicElements2D el;
    if (TwoBody::to_elements(mu, initial, el) != SolveStatus::Ok) {
        return false;
    }
    if (el.semi_latus_rectum <= 0.0 || !std::isfinite(el.semi_latus_rectum) ||
        el.eccentricity < 0.0 || !std::isfinite(el.eccentricity)) {
        return false;
    }

    const double e = el.eccentricity;
    const double r_min = el.semi_latus_rectum / (1.0 + e);
    double r_max = std::numeric_limits<double>::infinity();
    if (el.type == ConicType::Ellipse && e < 1.0) {
        const double denom = 1.0 - e;
        if (denom <= 0.0) return false;
        r_max = el.semi_latus_rectum / denom;
    }
    if (!std::isfinite(r_min) || r_min < 0.0) return false;

    out.min_radius = r_min;
    out.max_radius = r_max;
    return true;
}

bool child_soi_radially_overlaps(const RadialRange& range, const CachedChild& child,
                                 double root_eps) {
    const double child_inner = std::max(0.0, child.orbit_radius - child.soi_radius);
    const double child_outer = child.orbit_radius + child.soi_radius;
    return range.max_radius + root_eps >= child_inner &&
           range.min_radius - root_eps <= child_outer;
}

SolveStatus collect_active_children(const BodySystem& bodies, BodyId central_body,
                                    double mu, const State2& initial_state,
                                    double root_eps, ChildCache& out) {
    const BodyId* child_begin = bodies.children_begin(central_body);
    const BodyId* child_end = bodies.children_end(central_body);
    const size_t num_children = static_cast<size_t>(child_end - child_begin);
    if (num_children > kMaxChildrenPerBody) {
        return SolveStatus::InvalidInput;
    }

    out = ChildCache{};
    for (const BodyId* p = child_begin; p != child_end; ++p) {
        const BodyDef* ch = bodies.get_body(*p);
        if (ch == nullptr) continue;
        out.children[out.count++] =
            CachedChild{ch->id, ch->soi_radius, ch->orbit_radius,
                        ch->angular_rate, ch->phase_at_epoch};
    }

    if (conic_radial_range(mu, initial_state, out.radial_range)) {
        out.have_radial_range = true;
    }
    if (out.count > 0 && out.have_radial_range) {
        size_t active_count = 0;
        for (size_t i = 0; i < out.count; ++i) {
            if (child_soi_radially_overlaps(out.radial_range, out.children[i],
                                            root_eps)) {
                out.children[active_count++] = out.children[i];
            }
        }
        out.radially_culled = active_count == 0;
        out.count = active_count;
    }

    return SolveStatus::Ok;
}

void consider_start_boundary_events(const EventSearchRequest& req, BodyId parent_id,
                                    double body_radius, double soi_radius,
                                    const CachedChild* children, size_t child_count,
                                    double root_eps, double time_eps,
                                    CandidateSet& candidates) {
    const State2& initial = req.initial_state;

    // Impact at start: fire unconditionally if already at-or-inside the body radius.
    double f_impact = detail::impact_function(initial.r, body_radius);
    if (f_impact <= root_eps) {
        candidates.consider(make_event(EventType::Impact, req.start_time, req.central_body,
                                       InvalidBody, initial),
                            time_eps);
    }

    // SOI exit at start: fire only if outside by > root_eps OR on boundary with outward trend.
    double f_exit = detail::soi_exit_function(initial.r, soi_radius);
    double r_dot_v = dot(initial.r, initial.v);
    bool outside = f_exit > root_eps;
    bool on_boundary_outward = std::abs(f_exit) <= root_eps && r_dot_v > 0.0;
    if (outside || on_boundary_outward) {
        candidates.consider(make_event(EventType::SoiExit, req.start_time,
                                       req.central_body, parent_id, initial),
                            time_eps);
    }

    // Child SOI entry at start: fire only if inside by > root_eps OR on boundary with inward trend.
    for (size_t i = 0; i < child_count; ++i) {
        const CachedChild& child = children[i];
        Vec2 r_child = child_position(child, req.start_time);
        Vec2 v_child = child_velocity(child, req.start_time);
        Vec2 d = initial.r - r_child;
        Vec2 v_rel = initial.v - v_child;
        double f = detail::child_entry_function(initial.r, r_child, child.soi_radius);
        double d_dot_v = dot(d, v_rel);  // sign of d/dt(dist^2)/2 = sign of d/dt(dist)
        bool inside = f < -root_eps;
        bool on_boundary_inward = std::abs(f) <= root_eps && d_dot_v < 0.0;
        if (inside || on_boundary_inward) {
            candidates.consider(make_event(EventType::SoiEntry, req.start_time,
                                           req.central_body, child.id, initial),
                                time_eps);
        }
    }
}

SolveStatus consider_analytic_central_events(const EventSearchRequest& req,
                                             double mu, double body_radius,
                                             double soi_radius, BodyId parent_id,
                                             double horizon, double time_eps,
                                             const CachedConicPropagator& propagator,
                                             CandidateSet& candidates) {
    RadiusCrossing impact =
        first_central_radius_crossing(mu, req.initial_state, body_radius,
                                      RadiusCrossingDirection::Inward, horizon,
                                      time_eps);
    if (impact.found) {
        State2 s_event;
        SolveStatus ps = propagator.propagate(impact.dt, s_event);
        if (ps != SolveStatus::Ok) return ps;
        candidates.consider(make_event(EventType::Impact, req.start_time + impact.dt,
                                       req.central_body, InvalidBody, s_event),
                            time_eps);
    }

    RadiusCrossing exit =
        first_central_radius_crossing(mu, req.initial_state, soi_radius,
                                      RadiusCrossingDirection::Outward, horizon,
                                      time_eps);
    if (exit.found) {
        State2 s_event;
        SolveStatus ps = propagator.propagate(exit.dt, s_event);
        if (ps != SolveStatus::Ok) return ps;
        candidates.consider(make_event(EventType::SoiExit, req.start_time + exit.dt,
                                       req.central_body, parent_id, s_event),
                            time_eps);
    }

    return SolveStatus::Ok;
}

ScanConfig choose_scan_config(const EventSearchRequest& req, double horizon,
                              double body_radius, double soi_radius,
                              const CachedChild* children, size_t child_count,
                              bool have_radial_range, const RadialRange& radial_range,
                              double root_eps, double time_eps) {
    ScanConfig config;
    config.scan_impact =
        !have_radial_range || radial_range.min_radius <= body_radius + root_eps;
    config.scan_exit =
        !have_radial_range || radial_range.max_radius >= soi_radius - root_eps;

    config.spacecraft_speed = length(req.initial_state.v);
    config.max_relative_speed = config.spacecraft_speed;

    double min_child_soi = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < child_count; ++i) {
        const CachedChild& child = children[i];
        double v_moon = std::abs(child.orbit_radius * child.angular_rate);
        config.max_relative_speed =
            std::max(config.max_relative_speed, config.spacecraft_speed + v_moon);
        min_child_soi = std::min(min_child_soi, child.soi_radius);
    }

    const double v_floor = 0.1;
    double dt_child = std::isfinite(min_child_soi)
                          ? (min_child_soi / std::max(config.max_relative_speed, v_floor))
                          : std::numeric_limits<double>::infinity();
    double dt_default = horizon / 50.0;
    double dt_floor = std::max(time_eps * 100.0, 1e-12);
    config.dt = std::min(dt_child, dt_default);
    config.dt = std::max(config.dt, dt_floor);
    if (config.dt > horizon) config.dt = horizon;

    return config;
}

}  // namespace

EventDetector::EventDetector(const BodySystem& bodies, const Tolerances& tolerances,
                             EventDetectorDiagnostics* diagnostics)
    : bodies_(bodies), tolerances_(tolerances), diagnostics_(diagnostics) {}

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
    const CachedConicPropagator propagator(mu, req.initial_state);

    ChildCache child_cache;
    SolveStatus child_status =
        collect_active_children(bodies_, req.central_body, mu, req.initial_state,
                                root_eps, child_cache);
    if (child_status != SolveStatus::Ok) return child_status;
    const CachedChild* children = child_cache.children;
    const size_t child_count = child_cache.count;

    auto set_time_limit_output = [&]() {
        out = make_time_limit_event(req, propagator, horizon);
    };

    // --- Step 3: candidate arbitration ---
    CandidateSet candidates;

    // Root evaluation helpers. These operate on a pre-propagated spacecraft state.
    auto eval_impact = [&](const State2& s) { return length(s.r) - body_radius; };
    auto eval_exit = [&](const State2& s) { return length(s.r) - soi_radius; };
    auto eval_child = [&](const State2& s, size_t child_index, double abs_t) {
        const CachedChild& child = children[child_index];
        Vec2 rc = child_position(child, abs_t);
        return length(s.r - rc) - child.soi_radius;
    };

    // --- Step 4: start-state / boundary policy (spec 6.3) ---

    consider_start_boundary_events(req, parent_id, body_radius, soi_radius,
                                   children, child_count, root_eps, time_eps,
                                   candidates);

    // Degenerate horizon: skip scan entirely.
    if (horizon <= time_eps) {
        if (candidates.have_best) out = candidates.best;
        else set_time_limit_output();
        return SolveStatus::Ok;
    }

    // Central body impact and SOI exit depend only on the spacecraft's fixed conic
    // around the current central body, so solve their radius crossings analytically
    // before falling back to the coarse scan used for moving child SOIs.
    SolveStatus central_status =
        consider_analytic_central_events(req, mu, body_radius, soi_radius,
                                         parent_id, horizon, time_eps,
                                         propagator, candidates);
    if (central_status != SolveStatus::Ok) return central_status;

    if (child_cache.radially_culled) {
        if (candidates.have_best) out = candidates.best;
        else set_time_limit_output();
        return SolveStatus::Ok;
    }
    if (child_count == 0 && parent_id != InvalidBody && candidates.have_best) {
        out = candidates.best;
        return SolveStatus::Ok;
    }

    // --- Step 5: coarse scan step size ---
    // Adaptive step bounded by smallest child SOI / max relative speed, plus safety caps.
    ScanConfig scan = choose_scan_config(req, horizon, body_radius, soi_radius,
                                         children, child_count,
                                         child_cache.have_radial_range,
                                         child_cache.radial_range, root_eps,
                                         time_eps);

    // --- Step 6: coarse scan loop ---
    double t_lo = req.start_time;
    State2 s_lo = req.initial_state;
    double f_imp_lo = scan.scan_impact ? eval_impact(s_lo) : 0.0;
    double f_exit_lo = scan.scan_exit ? eval_exit(s_lo) : 0.0;
    double f_child_lo[kMaxChildrenPerBody];
    for (size_t i = 0; i < child_count; ++i) {
        f_child_lo[i] = eval_child(s_lo, i, t_lo);
    }

    // Hard cap on scan iterations prevents infinite loops from pathological inputs.
    const int max_scan_steps = 200000;
    int scan_steps = 0;

    while (t_lo < req.time_limit && scan_steps < max_scan_steps) {
        ++scan_steps;
        if (diagnostics_ != nullptr) ++diagnostics_->scan_steps;
        double t_hi = std::min(t_lo + scan.dt, req.time_limit);
        State2 s_hi;
        SolveStatus ps = propagator.propagate(t_hi - req.start_time, s_hi);
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
                SolveStatus ms = propagator.propagate(t_m - req.start_time, s_m);
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
                    (void)propagator.propagate(t - req.start_time, s);
                    return eval_impact(s);
                };
                double t_root = t_lo;
                if (diagnostics_ != nullptr) ++diagnostics_->root_refinements;
                SolveStatus rs = detail::refine_root_bisection(
                    f_of_t, t_lo, t_valid_lo, f_imp_lo, f_imp_end, root_eps, max_iter,
                    t_root);
                if (rs != SolveStatus::Ok) return rs;

                State2 s_root;
                (void)propagator.propagate(t_root - req.start_time, s_root);
                PredictedEvent ev;
                ev.type = EventType::Impact;
                ev.time = t_root;
                ev.from_body = req.central_body;
                ev.to_body = InvalidBody;
                ev.state = s_root;
                candidates.consider(ev, time_eps);
            }
            if (candidates.have_best) out = candidates.best;
            else set_time_limit_output();
            return SolveStatus::Ok;
        }

        double f_imp_hi = f_imp_lo;
        if (scan.scan_impact) {
            // --- Impact bracket ---
            f_imp_hi = eval_impact(s_hi);
            auto f_impact_of_t = [&](double t) {
                State2 s;
                (void)propagator.propagate(t - req.start_time, s);
                return eval_impact(s);
            };
            if (f_imp_lo > 0.0 && f_imp_hi <= 0.0) {
            double t_root = t_lo;
            if (diagnostics_ != nullptr) ++diagnostics_->root_refinements;
            SolveStatus rs = detail::refine_root_bisection(
                f_impact_of_t, t_lo, t_hi, f_imp_lo, f_imp_hi, root_eps, max_iter,
                t_root);
            if (rs != SolveStatus::Ok) return rs;

            State2 s_root;
            (void)propagator.propagate(t_root - req.start_time, s_root);
            PredictedEvent ev;
            ev.type = EventType::Impact;
            ev.time = t_root;
            ev.from_body = req.central_body;
            ev.to_body = InvalidBody;
            ev.state = s_root;
            candidates.consider(ev, time_eps);
        } else if (f_imp_lo > 0.0 && f_imp_hi > 0.0) {
            double t_mid = 0.5 * (t_lo + t_hi);
            double f_mid = f_impact_of_t(t_mid);
            double proximity =
                std::max({root_eps * 10.0, body_radius,
                          0.5 * scan.spacecraft_speed * scan.dt});
            bool possible_minimum = f_mid <= f_imp_lo && f_mid <= f_imp_hi;
            bool near_boundary = f_mid <= proximity;
            if (!(possible_minimum && near_boundary)) {
                // Same-sign windows are common. Only refine when the coarse
                // midpoint suggests a local minimum near the boundary.
            } else {
            double t_min = t_lo;
            double f_min = f_imp_lo;
            if (diagnostics_ != nullptr) ++diagnostics_->impact_minimum_refinements;
            SolveStatus ms = detail::refine_minimum(
                f_impact_of_t, t_lo, t_hi, root_eps, max_iter, t_min, f_min);
            if (ms != SolveStatus::Ok && f_min > root_eps) {
                // The bounded minimum search did not converge, but the best
                // sampled value is still comfortably outside the boundary.
                // Continue scanning instead of treating a non-candidate window
                // as an event-refinement failure.
                f_imp_lo = f_imp_hi;
            } else if (ms != SolveStatus::Ok) {
                return ms;
            }

            if (ms == SolveStatus::Ok && f_min <= root_eps) {
                double t_event = t_min;
                if (f_min <= 0.0 && t_min > t_lo + time_eps) {
                    double t_root = t_lo;
                    if (diagnostics_ != nullptr) ++diagnostics_->root_refinements;
                    SolveStatus rs = detail::refine_root_bisection(
                        f_impact_of_t, t_lo, t_min, f_imp_lo, f_min, root_eps, max_iter,
                        t_root);
                    if (rs != SolveStatus::Ok) return rs;
                    t_event = t_root;
                }

                State2 s_event;
                (void)propagator.propagate(t_event - req.start_time, s_event);
                PredictedEvent ev;
                ev.type = EventType::Impact;
                ev.time = t_event;
                ev.from_body = req.central_body;
                ev.to_body = InvalidBody;
                ev.state = s_event;
                candidates.consider(ev, time_eps);
            }
            }
        }
        }

        double f_exit_hi = f_exit_lo;
        if (scan.scan_exit) {
            // --- SOI exit bracket ---
            f_exit_hi = eval_exit(s_hi);
            auto f_exit_of_t = [&](double t) {
                State2 s;
                (void)propagator.propagate(t - req.start_time, s);
                return eval_exit(s);
            };
            if (f_exit_lo < 0.0 && f_exit_hi >= 0.0) {
            double t_root = t_lo;
            if (diagnostics_ != nullptr) ++diagnostics_->root_refinements;
            SolveStatus rs = detail::refine_root_bisection(
                f_exit_of_t, t_lo, t_hi, f_exit_lo, f_exit_hi, root_eps, max_iter,
                t_root);
            if (rs != SolveStatus::Ok) return rs;

            State2 s_root;
            (void)propagator.propagate(t_root - req.start_time, s_root);
            PredictedEvent ev;
            ev.type = EventType::SoiExit;
            ev.time = t_root;
            ev.from_body = req.central_body;
            ev.to_body = parent_id;
            ev.state = s_root;
            candidates.consider(ev, time_eps);
        } else if (f_exit_lo < 0.0 && f_exit_hi < 0.0) {
            double t_mid = 0.5 * (t_lo + t_hi);
            double f_mid = f_exit_of_t(t_mid);
            double proximity = std::max(root_eps * 10.0, soi_radius * 0.05);
            bool possible_maximum = f_mid >= f_exit_lo && f_mid >= f_exit_hi;
            bool near_boundary = f_mid >= -proximity;
            if (!(possible_maximum && near_boundary)) {
                // Same-sign windows are common. Only refine when the coarse
                // midpoint suggests a local maximum near the SOI boundary.
            } else {
            auto neg_exit_of_t = [&](double t) { return -f_exit_of_t(t); };
            double t_max = t_lo;
            double neg_f_max = -f_exit_lo;
            if (diagnostics_ != nullptr) ++diagnostics_->soi_exit_minimum_refinements;
            SolveStatus ms = detail::refine_minimum(
                neg_exit_of_t, t_lo, t_hi, root_eps, max_iter, t_max, neg_f_max);

            double f_max = -neg_f_max;
            if (ms != SolveStatus::Ok && f_max < -root_eps) {
                // As above, non-convergence in a non-candidate window should not
                // mask ordinary monotonic scan progress.
                f_exit_lo = f_exit_hi;
            } else if (ms != SolveStatus::Ok) {
                return ms;
            }

            if (ms == SolveStatus::Ok && f_max >= -root_eps) {
                double t_event = t_max;
                if (f_max >= 0.0 && t_max > t_lo + time_eps) {
                    double t_root = t_lo;
                    if (diagnostics_ != nullptr) ++diagnostics_->root_refinements;
                    SolveStatus rs = detail::refine_root_bisection(
                        f_exit_of_t, t_lo, t_max, f_exit_lo, f_max, root_eps, max_iter,
                        t_root);
                    if (rs != SolveStatus::Ok) return rs;
                    t_event = t_root;
                }

                State2 s_event;
                (void)propagator.propagate(t_event - req.start_time, s_event);
                PredictedEvent ev;
                ev.type = EventType::SoiExit;
                ev.time = t_event;
                ev.from_body = req.central_body;
                ev.to_body = parent_id;
                ev.state = s_event;
                candidates.consider(ev, time_eps);
            }
            }
        }
        }

        // --- Child entry brackets ---
        for (size_t ci = 0; ci < child_count; ++ci) {
            const CachedChild& child = children[ci];
            double f_ch_lo = f_child_lo[ci];
            double f_ch_hi = eval_child(s_hi, ci, t_hi);
            auto f_of_t = [&, ci](double t) {
                State2 s;
                (void)propagator.propagate(t - req.start_time, s);
                return eval_child(s, ci, t);
            };
            if (f_ch_lo > 0.0 && f_ch_hi <= 0.0) {
                double t_root = t_lo;
                if (diagnostics_ != nullptr) ++diagnostics_->root_refinements;
                SolveStatus rs = detail::refine_root_bisection(
                    f_of_t, t_lo, t_hi, f_ch_lo, f_ch_hi, root_eps, max_iter, t_root);
                if (rs != SolveStatus::Ok) return rs;

                State2 s_root;
                (void)propagator.propagate(t_root - req.start_time, s_root);
                PredictedEvent ev;
                ev.type = EventType::SoiEntry;
                ev.time = t_root;
                ev.from_body = req.central_body;
                ev.to_body = child.id;
                ev.state = s_root;
                candidates.consider(ev, time_eps);
            } else if (f_ch_lo > 0.0 && f_ch_hi > 0.0) {
                const double half_window_motion =
                    0.5 * scan.max_relative_speed * (t_hi - t_lo);
                double proximity =
                    std::max({root_eps * 10.0, child.soi_radius * 0.05,
                              half_window_motion});
                bool too_far_for_midpoint =
                    std::min(f_ch_lo, f_ch_hi) > proximity + half_window_motion;
                if (too_far_for_midpoint) {
                    // The distance function cannot plausibly dip close enough
                    // to the child SOI within this window to warrant a midpoint
                    // propagation.
                } else {
                double t_mid = 0.5 * (t_lo + t_hi);
                double f_mid = f_of_t(t_mid);
                bool possible_minimum = f_mid <= f_ch_lo && f_mid <= f_ch_hi;
                bool near_boundary = f_mid <= proximity;
                if (!(possible_minimum && near_boundary)) {
                    // Same-sign windows are common. Only refine when the coarse
                    // midpoint suggests a local minimum near the child SOI.
                } else {
                double t_min = t_lo;
                double f_min = f_ch_lo;
                if (diagnostics_ != nullptr) ++diagnostics_->child_entry_minimum_refinements;
                SolveStatus ms = detail::refine_minimum(
                    f_of_t, t_lo, t_hi, root_eps, max_iter, t_min, f_min);
                if (ms != SolveStatus::Ok) return ms;

                if (ms == SolveStatus::Ok && f_min <= root_eps) {
                    double t_event = t_min;
                    if (f_min <= 0.0 && t_min > t_lo + time_eps) {
                        double t_root = t_lo;
                        if (diagnostics_ != nullptr) ++diagnostics_->root_refinements;
                        SolveStatus rs = detail::refine_root_bisection(
                            f_of_t, t_lo, t_min, f_ch_lo, f_min, root_eps, max_iter,
                            t_root);
                        if (rs == SolveStatus::Ok) {
                            t_event = t_root;
                        } else {
                            return rs;
                        }
                    }
                    State2 s_event;
                    (void)propagator.propagate(t_event - req.start_time, s_event);
                    PredictedEvent ev;
                    ev.type = EventType::SoiEntry;
                    ev.time = t_event;
                    ev.from_body = req.central_body;
                    ev.to_body = child.id;
                    ev.state = s_event;
                    candidates.consider(ev, time_eps);
                }
                }
                }
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
        if (candidates.have_best && candidates.best.time < t_lo - time_eps) {
            break;
        }
    }

    if (candidates.have_best) {
        out = candidates.best;
    } else {
        set_time_limit_output();
    }
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
