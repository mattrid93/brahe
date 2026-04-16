#include "brahe/two_body.h"

#include "brahe/vec2.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace brahe {

static constexpr Tolerances kTol{};
static constexpr double kCircularThreshold = 1e-10;
static constexpr double kKeplerTol = 1e-14;

// ============================================================
// detail namespace — numerical helpers and anomaly conversions
// ============================================================

namespace detail {

double wrap_angle(double a) {
    a = std::fmod(a, 2.0 * M_PI);
    if (a > M_PI) a -= 2.0 * M_PI;
    if (a <= -M_PI) a += 2.0 * M_PI;
    return a;
}

double safe_acos(double x) {
    return std::acos(std::clamp(x, -1.0, 1.0));
}

double safe_asin(double x) {
    return std::asin(std::clamp(x, -1.0, 1.0));
}

// --- Elliptic anomaly conversions (0 <= e < 1) ---

double true_to_eccentric_anomaly(double nu, double e) {
    double sin_nu = std::sin(nu);
    double cos_nu = std::cos(nu);
    double sqrt_1me2 = std::sqrt(std::max(0.0, 1.0 - e * e));
    return std::atan2(sqrt_1me2 * sin_nu, e + cos_nu);
}

double eccentric_to_true_anomaly(double E, double e) {
    double sin_E = std::sin(E);
    double cos_E = std::cos(E);
    double sqrt_1me2 = std::sqrt(std::max(0.0, 1.0 - e * e));
    return std::atan2(sqrt_1me2 * sin_E, cos_E - e);
}

double eccentric_to_mean_anomaly(double E, double e) {
    return E - e * std::sin(E);
}

SolveStatus mean_to_eccentric_anomaly(double M, double e, double& out_E, int max_iterations) {
    // Newton's method: f(E) = E - e*sin(E) - M, f'(E) = 1 - e*cos(E)
    double E = M + e * std::sin(M); // initial guess (first-order correction)
    double tol = kKeplerTol * std::max(1.0, std::abs(M));
    for (int i = 0; i < max_iterations; ++i) {
        double f = E - e * std::sin(E) - M;
        if (std::abs(f) < tol) {
            out_E = E;
            return SolveStatus::Ok;
        }
        double fp = 1.0 - e * std::cos(E);
        if (std::abs(fp) < 1e-30) return SolveStatus::NumericalFailure;
        E -= f / fp;
    }
    return SolveStatus::NoConvergence;
}

// --- Hyperbolic anomaly conversions (e > 1) ---

double true_to_hyperbolic_anomaly(double nu, double e) {
    // tanh(H/2) = sqrt((e-1)/(e+1)) * tan(nu/2)
    // Use atan2-based formula for robustness:
    double half_nu = nu / 2.0;
    double sin_half = std::sin(half_nu);
    double cos_half = std::cos(half_nu);
    double factor = std::sqrt(std::max(0.0, (e - 1.0) / (e + 1.0)));
    double tanh_half_H = factor * sin_half / cos_half;

    // H = 2 * atanh(tanh_half_H)
    // atanh(x) = 0.5 * log((1+x)/(1-x))
    tanh_half_H = std::clamp(tanh_half_H, -1.0 + 1e-15, 1.0 - 1e-15);
    return 2.0 * std::atanh(tanh_half_H);
}

double hyperbolic_to_true_anomaly(double H, double e) {
    double half_H = H / 2.0;
    double sqrt_ep1 = std::sqrt(e + 1.0);
    double sqrt_em1 = std::sqrt(std::max(0.0, e - 1.0));
    return 2.0 * std::atan2(sqrt_ep1 * std::sinh(half_H), sqrt_em1 * std::cosh(half_H));
}

double hyperbolic_to_mean_anomaly(double H, double e) {
    return e * std::sinh(H) - H;
}

SolveStatus mean_to_hyperbolic_anomaly(double Mh, double e, double& out_H, int max_iterations) {
    // Newton's method: f(H) = e*sinh(H) - H - Mh, f'(H) = e*cosh(H) - 1
    // Initial guess
    double H;
    if (std::abs(Mh) < 1.0) {
        H = Mh / (e - 1.0);
    } else {
        H = std::copysign(std::log(2.0 * std::abs(Mh) / e + 1.0), Mh);
    }

    double tol = kKeplerTol * std::max(1.0, std::abs(Mh));
    for (int i = 0; i < max_iterations; ++i) {
        double f = e * std::sinh(H) - H - Mh;
        if (std::abs(f) < tol) {
            out_H = H;
            return SolveStatus::Ok;
        }
        double fp = e * std::cosh(H) - 1.0;
        if (std::abs(fp) < 1e-30) return SolveStatus::NumericalFailure;
        H -= f / fp;
    }
    return SolveStatus::NoConvergence;
}

// --- Radial (h = 0) Kepler solvers ---
//
// Purely radial trajectories degenerate the conic parameterization (p = 0)
// and must be propagated directly from the 1-D radial Kepler problem.
// Parametric forms used below:
//   Bound  (eps < 0): r = a*(1 - cos(u)),  t - t_peri = sqrt(a^3/mu)*(u - sin(u)), a = -mu/(2 eps) > 0
//   Hyper  (eps > 0): r = a*(cosh(F) - 1), t - t_peri = sqrt(a^3/mu)*(sinh(F) - F), a =  mu/(2 eps) > 0
//   Parab  (eps = 0): r^{3/2} = r_0^{3/2} + sign(v_r) * (3/2) sqrt(2 mu) * dt
// In the bound and hyperbolic forms the point r = 0 is the central-body singularity;
// propagation that would cross it fails with NumericalFailure.

SolveStatus solve_radial_elliptic_kepler(double M_reduced, double& out_u,
                                         int max_iterations = 100) {
    // Solve u - sin(u) = M_reduced for u in [0, 2*pi), with M_reduced in [0, 2*pi).
    double u;
    if (M_reduced < 1e-3) {
        u = std::cbrt(6.0 * M_reduced);                    // cube-root expansion near u=0
    } else if (M_reduced > 2.0 * M_PI - 1e-3) {
        u = 2.0 * M_PI - std::cbrt(6.0 * (2.0 * M_PI - M_reduced));
    } else {
        u = M_reduced + 0.5 * std::sin(M_reduced);
    }
    double tol = kKeplerTol * std::max(1.0, std::abs(M_reduced));
    for (int i = 0; i < max_iterations; ++i) {
        double f = u - std::sin(u) - M_reduced;
        if (std::abs(f) < tol) {
            out_u = u;
            return SolveStatus::Ok;
        }
        double fp = 1.0 - std::cos(u);
        if (fp < 1e-15) {
            // Near u = 0 or 2*pi (central singularity); a small nudge in the
            // direction of the residual keeps Newton off the flat spot.
            u += (f < 0.0 ? 1e-3 : -1e-3);
            continue;
        }
        double du = f / fp;
        if (du > 0.5) du = 0.5;
        if (du < -0.5) du = -0.5;
        u -= du;
        if (u < 0.0) u = 0.0;
        if (u > 2.0 * M_PI) u = 2.0 * M_PI;
    }
    return SolveStatus::NoConvergence;
}

SolveStatus solve_radial_hyperbolic_kepler(double M, double& out_F,
                                           int max_iterations = 100) {
    // Solve sinh(F) - F = M for F.  sign(F) = sign(M).
    double F;
    if (std::abs(M) < 1.0) {
        F = std::cbrt(6.0 * M);                            // small-M expansion
    } else {
        F = std::copysign(std::log(2.0 * std::abs(M) + 1.0), M);
    }
    double tol = kKeplerTol * std::max(1.0, std::abs(M));
    for (int i = 0; i < max_iterations; ++i) {
        double f = std::sinh(F) - F - M;
        if (std::abs(f) < tol) {
            out_F = F;
            return SolveStatus::Ok;
        }
        double fp = std::cosh(F) - 1.0;
        if (fp < 1e-15) {
            F += (f < 0.0 ? 1e-3 : -1e-3);
            continue;
        }
        F -= f / fp;
    }
    return SolveStatus::NoConvergence;
}

SolveStatus radial_propagate(double mu, const State2& initial, double dt,
                             State2& out_state) {
    double r_0 = length(initial.r);
    Vec2 u_hat = {initial.r.x / r_0, initial.r.y / r_0};
    double v_r_0 = initial.v.x * u_hat.x + initial.v.y * u_hat.y;   // signed radial speed
    double eps = 0.5 * v_r_0 * v_r_0 - mu / r_0;
    // Band width (parabolic) scaled by characteristic energy mu/r_0.
    double eps_band = kTol.parabolic_eccentricity_band * (mu / r_0);

    double r_new = 0.0;
    double v_r_new = 0.0;

    if (std::abs(eps) <= eps_band) {
        // --- Parabolic radial ---
        // r^{3/2} = r_0^{3/2} + sign(v_r_0) * 1.5 * sqrt(2 mu) * dt
        double sign_v = (v_r_0 >= 0.0) ? 1.0 : -1.0;
        double r_3_2 = std::pow(r_0, 1.5) + sign_v * 1.5 * std::sqrt(2.0 * mu) * dt;
        if (r_3_2 <= 0.0) {
            return SolveStatus::NumericalFailure;         // inbound crossed r=0
        }
        r_new = std::pow(r_3_2, 2.0 / 3.0);
        v_r_new = sign_v * std::sqrt(2.0 * mu / r_new);
    } else if (eps < 0.0) {
        // --- Bound radial ---
        double a = -mu / (2.0 * eps);                     // a > 0
        double cos_E0 = std::clamp(1.0 - r_0 / a, -1.0, 1.0);
        double E_0 = std::acos(cos_E0);                   // [0, pi]
        if (v_r_0 < 0.0) E_0 = 2.0 * M_PI - E_0;           // inbound: E_0 in (pi, 2pi)
        double M_0 = E_0 - std::sin(E_0);
        double n = std::sqrt(mu / (a * a * a));
        double M_new = M_0 + n * dt;
        if (M_new < 0.0 || M_new >= 2.0 * M_PI) {
            // Propagation would cross the r = 0 central-body singularity.
            return SolveStatus::NumericalFailure;
        }
        double E_new;
        SolveStatus ks = solve_radial_elliptic_kepler(M_new, E_new);
        if (ks != SolveStatus::Ok) return ks;
        double one_minus_cos = 1.0 - std::cos(E_new);
        r_new = a * one_minus_cos;
        if (one_minus_cos < 1e-30) return SolveStatus::NumericalFailure;
        v_r_new = std::sqrt(mu / a) * std::sin(E_new) / one_minus_cos;
    } else {
        // --- Hyperbolic radial ---
        double a = mu / (2.0 * eps);                      // a > 0
        double cosh_F0 = 1.0 + r_0 / a;
        double F_0 = std::acosh(std::max(1.0, cosh_F0));  // >= 0
        if (v_r_0 < 0.0) F_0 = -F_0;
        double M_0 = std::sinh(F_0) - F_0;
        double n = std::sqrt(mu / (a * a * a));
        double M_new = M_0 + n * dt;
        if (F_0 < 0.0 && M_new >= 0.0) {
            // Inbound trajectory reaches r = 0 within dt.
            return SolveStatus::NumericalFailure;
        }
        double F_new;
        SolveStatus ks = solve_radial_hyperbolic_kepler(M_new, F_new);
        if (ks != SolveStatus::Ok) return ks;
        double cosh_m1 = std::cosh(F_new) - 1.0;
        if (cosh_m1 < 1e-30) return SolveStatus::NumericalFailure;
        r_new = a * cosh_m1;
        v_r_new = std::sqrt(mu / a) * std::sinh(F_new) / cosh_m1;
    }

    out_state.r = {r_new * u_hat.x, r_new * u_hat.y};
    out_state.v = {v_r_new * u_hat.x, v_r_new * u_hat.y};
    return SolveStatus::Ok;
}

} // namespace detail

// ============================================================
// Input validation helper
// ============================================================

static SolveStatus validate_mu_state(double mu, const State2& state) {
    if (mu <= 0.0 || !std::isfinite(mu)) return SolveStatus::InvalidInput;
    if (!is_finite(state.r) || !is_finite(state.v)) return SolveStatus::InvalidInput;
    if (length(state.r) < kTol.position_epsilon) return SolveStatus::InvalidInput;
    return SolveStatus::Ok;
}

// ============================================================
// TwoBody — orbital invariants
// ============================================================

double TwoBody::specific_energy(double mu, const State2& state) {
    double v_sq = length_squared(state.v);
    double r_mag = length(state.r);
    return v_sq / 2.0 - mu / r_mag;
}

double TwoBody::specific_angular_momentum_z(const State2& state) {
    return state.r.x * state.v.y - state.r.y * state.v.x;
}

double TwoBody::eccentricity(double mu, const State2& state) {
    double r_mag = length(state.r);
    double v_sq = length_squared(state.v);
    double rdotv = dot(state.r, state.v);

    // e_vec = (1/mu) * [(v² - mu/r)*r - dot(r,v)*v]
    Vec2 e_vec = (1.0 / mu) * ((v_sq - mu / r_mag) * state.r - rdotv * state.v);
    return length(e_vec);
}

double TwoBody::semi_major_axis(double mu, const State2& state) {
    double eps = specific_energy(mu, state);
    if (std::abs(eps) < 1e-14) {
        return std::numeric_limits<double>::infinity();
    }
    return -mu / (2.0 * eps);
}

double TwoBody::semi_latus_rectum(double mu, const State2& state) {
    double h = specific_angular_momentum_z(state);
    return h * h / mu;
}

double TwoBody::periapsis_radius(double mu, const State2& state) {
    double p = semi_latus_rectum(mu, state);
    double e = eccentricity(mu, state);
    return p / (1.0 + e);
}

double TwoBody::apoapsis_radius(double mu, const State2& state) {
    double e = eccentricity(mu, state);
    if (e >= 1.0 - kTol.parabolic_eccentricity_band) {
        return std::numeric_limits<double>::infinity();
    }
    double p = semi_latus_rectum(mu, state);
    return p / (1.0 - e);
}

// ============================================================
// TwoBody — classification
// ============================================================

SolveStatus TwoBody::classify(double mu, const State2& state, ConicType& out_type) {
    SolveStatus s = validate_mu_state(mu, state);
    if (s != SolveStatus::Ok) return s;

    double e = eccentricity(mu, state);
    double band = kTol.parabolic_eccentricity_band;

    if (e < 1.0 - band) {
        out_type = ConicType::Ellipse;
    } else if (e > 1.0 + band) {
        out_type = ConicType::Hyperbola;
    } else {
        out_type = ConicType::ParabolaLike;
    }
    return SolveStatus::Ok;
}

// ============================================================
// TwoBody — Cartesian -> conic elements
// ============================================================

SolveStatus TwoBody::to_elements(double mu, const State2& state, ConicElements2D& out) {
    SolveStatus s = validate_mu_state(mu, state);
    if (s != SolveStatus::Ok) return s;

    double r_mag = length(state.r);
    double v_sq = length_squared(state.v);
    double rdotv = dot(state.r, state.v);
    double h = specific_angular_momentum_z(state);
    double eps = v_sq / 2.0 - mu / r_mag;
    double p = h * h / mu;

    // Eccentricity vector
    Vec2 e_vec = (1.0 / mu) * ((v_sq - mu / r_mag) * state.r - rdotv * state.v);
    double e = length(e_vec);

    // Classification
    double band = kTol.parabolic_eccentricity_band;
    ConicType type;
    if (e < 1.0 - band) {
        type = ConicType::Ellipse;
    } else if (e > 1.0 + band) {
        type = ConicType::Hyperbola;
    } else {
        type = ConicType::ParabolaLike;
    }

    // Semi-major axis
    double a;
    if (std::abs(eps) < 1e-14) {
        a = std::numeric_limits<double>::infinity();
    } else {
        a = -mu / (2.0 * eps);
    }

    // Argument of periapsis and true anomaly
    double omega, nu;
    if (e < kCircularThreshold) {
        // Circular convention: omega = 0, nu = position angle
        omega = 0.0;
        nu = std::atan2(state.r.y, state.r.x);
    } else {
        omega = std::atan2(e_vec.y, e_vec.x);
        double theta = std::atan2(state.r.y, state.r.x);
        nu = detail::wrap_angle(theta - omega);
    }

    // Fill output
    out.type = type;
    out.mu = mu;
    out.semi_major_axis = a;
    out.eccentricity = e;
    out.semi_latus_rectum = p;
    out.angular_momentum_z = h;
    out.argument_of_periapsis = detail::wrap_angle(omega);
    out.true_anomaly = detail::wrap_angle(nu);

    // Anomaly fields
    out.mean_anomaly = 0.0;
    out.eccentric_anomaly = 0.0;
    out.hyperbolic_anomaly = 0.0;

    if (type == ConicType::Ellipse) {
        double E = detail::true_to_eccentric_anomaly(nu, e);
        double M = detail::eccentric_to_mean_anomaly(E, e);
        out.eccentric_anomaly = E;
        out.mean_anomaly = M;
    } else if (type == ConicType::Hyperbola) {
        double H = detail::true_to_hyperbolic_anomaly(nu, e);
        double Mh = detail::hyperbolic_to_mean_anomaly(H, e);
        out.hyperbolic_anomaly = H;
        out.mean_anomaly = Mh;
    } else {
        // ParabolaLike: store parabolic mean anomaly D + D^3/3
        double D = std::tan(nu / 2.0);
        out.mean_anomaly = D + D * D * D / 3.0;
    }

    return SolveStatus::Ok;
}

// ============================================================
// TwoBody — conic elements -> Cartesian
// ============================================================

SolveStatus TwoBody::from_elements(const ConicElements2D& el, State2& out_state) {
    if (el.mu <= 0.0 || !std::isfinite(el.mu)) return SolveStatus::InvalidInput;
    if (el.semi_latus_rectum <= 0.0 || !std::isfinite(el.semi_latus_rectum))
        return SolveStatus::InvalidInput;
    if (el.eccentricity < 0.0 || !std::isfinite(el.eccentricity))
        return SolveStatus::InvalidInput;

    double mu = el.mu;
    double p = el.semi_latus_rectum;
    double e = el.eccentricity;
    double nu = el.true_anomaly;
    double omega = el.argument_of_periapsis;

    double cos_nu = std::cos(nu);
    double sin_nu = std::sin(nu);
    double denom = 1.0 + e * cos_nu;
    if (std::abs(denom) < 1e-15) return SolveStatus::NumericalFailure;

    double r_mag = p / denom;

    // Position in perifocal frame
    Vec2 r_pf = {r_mag * cos_nu, r_mag * sin_nu};

    // Velocity in perifocal frame
    // v_pf = (mu/h) * [-sin(nu), e + cos(nu)]
    // where h = angular_momentum_z (signed), |h| = sqrt(mu * p)
    double abs_h = std::sqrt(mu * p);
    double h = (el.angular_momentum_z < 0.0) ? -abs_h : abs_h;
    if (std::abs(h) < 1e-30) return SolveStatus::NumericalFailure;

    Vec2 v_pf = {(mu / h) * (-sin_nu), (mu / h) * (e + cos_nu)};

    // Rotate from perifocal to inertial by argument_of_periapsis
    double cos_w = std::cos(omega);
    double sin_w = std::sin(omega);

    out_state.r.x = cos_w * r_pf.x - sin_w * r_pf.y;
    out_state.r.y = sin_w * r_pf.x + cos_w * r_pf.y;
    out_state.v.x = cos_w * v_pf.x - sin_w * v_pf.y;
    out_state.v.y = sin_w * v_pf.x + cos_w * v_pf.y;

    return SolveStatus::Ok;
}

// ============================================================
// TwoBody — analytic propagation
// ============================================================

SolveStatus TwoBody::propagate(double mu, const State2& initial, double dt,
                                State2& out_state) {
    // Input validation
    if (mu <= 0.0 || !std::isfinite(mu)) return SolveStatus::InvalidInput;
    if (!is_finite(initial.r) || !is_finite(initial.v)) return SolveStatus::InvalidInput;
    if (length(initial.r) < kTol.position_epsilon) return SolveStatus::InvalidInput;
    if (dt < 0.0 || !std::isfinite(dt)) return SolveStatus::InvalidInput;

    if (dt == 0.0) {
        out_state = initial;
        return SolveStatus::Ok;
    }

    // Radial (h = 0) dispatch.  The conic parameterization r = p / (1 + e cos nu)
    // degenerates at p = 0, so we propagate the 1-D radial Kepler problem directly.
    // Detection uses a relative tolerance scaled by the problem's characteristic
    // value of |r x v| to avoid false positives on ordinary low-h orbits.
    {
        double r_mag = length(initial.r);
        double v_mag = length(initial.v);
        double h = specific_angular_momentum_z(initial);
        double h_threshold = 1e-12 * std::max(1.0, r_mag * v_mag);
        if (std::abs(h) <= h_threshold) {
            return detail::radial_propagate(mu, initial, dt, out_state);
        }
    }

    // Convert to elements
    ConicElements2D el;
    SolveStatus s = to_elements(mu, initial, el);
    if (s != SolveStatus::Ok) return s;

    double e = el.eccentricity;
    double p = el.semi_latus_rectum;
    double h = el.angular_momentum_z;
    double sign_h = (h >= 0.0) ? 1.0 : -1.0;
    double band = kTol.parabolic_eccentricity_band;

    double nu_new;

    if (std::abs(e - 1.0) <= band) {
        // --- Parabolic path (Barker's equation) ---
        double D0 = std::tan(el.true_anomaly / 2.0);
        double coeff = std::sqrt(2.0 * mu / (p * p * p));
        double W = D0 + D0 * D0 * D0 / 3.0 + sign_h * coeff * dt;

        // Solve D + D^3/3 = W analytically
        // Depressed cubic: D^3 + 3D - 3W = 0
        // Solution: D = 2*sinh(asinh(3W/2)/3)
        double z = 1.5 * W;
        double D1 = 2.0 * std::sinh(std::asinh(z) / 3.0);
        nu_new = 2.0 * std::atan(D1);
    } else if (e < 1.0) {
        // --- Elliptic path ---
        double a = el.semi_major_axis;
        if (a <= 0.0 || !std::isfinite(a)) return SolveStatus::NumericalFailure;
        double n = std::sqrt(mu / (a * a * a));
        double M_new = el.mean_anomaly + sign_h * n * dt;

        double E_new;
        SolveStatus ks = detail::mean_to_eccentric_anomaly(M_new, e, E_new);
        if (ks != SolveStatus::Ok) return ks;
        nu_new = detail::eccentric_to_true_anomaly(E_new, e);
    } else {
        // --- Hyperbolic path ---
        double a = el.semi_major_axis; // negative for hyperbola
        if (a >= 0.0 || !std::isfinite(a)) return SolveStatus::NumericalFailure;
        double n = std::sqrt(mu / (-a * a * a));
        double Mh_new = el.mean_anomaly + sign_h * n * dt;

        double H_new;
        SolveStatus ks = detail::mean_to_hyperbolic_anomaly(Mh_new, e, H_new);
        if (ks != SolveStatus::Ok) return ks;
        nu_new = detail::hyperbolic_to_true_anomaly(H_new, e);
    }

    // Update elements and reconstruct Cartesian state
    el.true_anomaly = nu_new;
    return from_elements(el, out_state);
}

} // namespace brahe
