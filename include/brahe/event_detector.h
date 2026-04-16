#pragma once

#include "brahe/body_system.h"
#include "brahe/types.h"
#include "brahe/vec2.h"

namespace brahe {

// --- Public event types ---

struct PredictedEvent {
    EventType type = EventType::None;
    double time = 0.0;
    BodyId from_body = InvalidBody;
    BodyId to_body = InvalidBody;  // child for SoiEntry, parent for SoiExit, InvalidBody otherwise
    State2 state;                  // spacecraft state in central-body frame at event time (pre-event)
};

struct EventSearchRequest {
    BodyId central_body = InvalidBody;
    double start_time = 0.0;
    State2 initial_state;  // state relative to central_body at start_time
    double time_limit = 0.0;
};

struct EventDetectorDiagnostics {
    int scan_steps = 0;
    int root_refinements = 0;
    int impact_minimum_refinements = 0;
    int soi_exit_minimum_refinements = 0;
    int child_entry_minimum_refinements = 0;
};

// --- Public detector API ---

class EventDetector {
public:
    EventDetector(const BodySystem& bodies, const Tolerances& tolerances = {},
                  EventDetectorDiagnostics* diagnostics = nullptr);

    SolveStatus find_next_event(const EventSearchRequest& req,
                                PredictedEvent& out_event) const;

private:
    const BodySystem& bodies_;
    Tolerances tolerances_;
    EventDetectorDiagnostics* diagnostics_ = nullptr;
};

// --- Internal helpers exposed for testing ---

namespace detail {

// Root functions: distance to boundary.
// Positive = outside boundary, negative = inside, zero = on boundary.

inline double impact_function(const Vec2& r, double body_radius) {
    return length(r) - body_radius;
}

inline double soi_exit_function(const Vec2& r, double soi_radius) {
    return length(r) - soi_radius;
}

inline double child_entry_function(const Vec2& spacecraft_pos, const Vec2& child_pos,
                                   double child_soi_radius) {
    return length(spacecraft_pos - child_pos) - child_soi_radius;
}

// Bounded bisection root refinement.
// Requires: f_lo and f_hi have opposite signs.
// F signature: double f(double t)
// Returns refined root time in out_t.
template <typename F>
SolveStatus refine_root_bisection(F&& eval, double t_lo, double t_hi, double f_lo, double f_hi,
                                  double tol, int max_iterations, double& out_t) {
    if (max_iterations <= 0) {
        out_t = 0.5 * (t_lo + t_hi);
        return SolveStatus::NoConvergence;
    }
    if (t_lo >= t_hi) {
        out_t = 0.5 * (t_lo + t_hi);
        return SolveStatus::NoConvergence;
    }

    // Determine the sign convention at the two endpoints so we can narrow the
    // bracket correctly without swapping t_lo / t_hi (which would break the
    // t_lo < t_hi invariant that the convergence check depends on).
    //
    // polarity = +1 when f_lo <= 0 <= f_hi (f increasing through the root)
    // polarity = -1 when f_lo >= 0 >= f_hi (f decreasing through the root)
    double polarity;
    if (f_lo <= 0.0 && f_hi >= 0.0) {
        polarity = 1.0;
    } else if (f_lo >= 0.0 && f_hi <= 0.0) {
        polarity = -1.0;
    } else {
        // Endpoints do not bracket a sign change.
        out_t = 0.5 * (t_lo + t_hi);
        return SolveStatus::NoConvergence;
    }

    for (int i = 0; i < max_iterations; ++i) {
        double t_mid = 0.5 * (t_lo + t_hi);
        if (t_hi - t_lo <= tol) {
            out_t = t_mid;
            return SolveStatus::Ok;
        }
        double f_mid = eval(t_mid);
        // Keep t_lo as the "low-polarity-signed" end. When polarity == +1 that
        // means f_lo is on the negative side; when polarity == -1 it means f_lo
        // is on the positive side. Either way, the root is in [t_mid, t_hi]
        // iff polarity * f_mid <= 0.
        if (polarity * f_mid <= 0.0) {
            t_lo = t_mid;
            f_lo = f_mid;
        } else {
            t_hi = t_mid;
            f_hi = f_mid;
        }
    }
    out_t = 0.5 * (t_lo + t_hi);
    return SolveStatus::NoConvergence;
}

// Bounded golden-section minimum search.
// Finds minimum of f in [t_lo, t_hi].
// F signature: double f(double t)
// Returns minimum time in out_t, minimum value in out_f_min.
template <typename F>
SolveStatus refine_minimum(F&& eval, double t_lo, double t_hi, double tol, int max_iterations,
                           double& out_t, double& out_f_min) {
    if (max_iterations <= 0) {
        out_t = 0.5 * (t_lo + t_hi);
        out_f_min = eval(out_t);
        return SolveStatus::NoConvergence;
    }

    constexpr double phi = 0.6180339887498949;  // (sqrt(5) - 1) / 2

    double a = t_lo;
    double b = t_hi;
    double c = b - phi * (b - a);
    double d = a + phi * (b - a);
    double fc = eval(c);
    double fd = eval(d);

    for (int i = 0; i < max_iterations; ++i) {
        if (b - a <= tol) {
            break;
        }
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - phi * (b - a);
            fc = eval(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + phi * (b - a);
            fd = eval(d);
        }
    }

    if (b - a > tol) {
        out_t = 0.5 * (a + b);
        out_f_min = eval(out_t);
        return SolveStatus::NoConvergence;
    }

    out_t = 0.5 * (a + b);
    out_f_min = eval(out_t);
    return SolveStatus::Ok;
}

// Deterministic event comparison.
// Returns true if event a should be selected over event b.
bool event_precedes(const PredictedEvent& a, const PredictedEvent& b, double time_epsilon);

}  // namespace detail

}  // namespace brahe
