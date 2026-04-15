#pragma once

#include "brahe/conics.h"
#include "brahe/types.h"

namespace brahe {

class TwoBody {
public:
    // Classification
    static SolveStatus classify(double mu, const State2& state, ConicType& out_type);

    // Cartesian <-> conic element conversion
    static SolveStatus to_elements(double mu, const State2& state, ConicElements2D& out);
    static SolveStatus from_elements(const ConicElements2D& elements, State2& out_state);

    // Analytic propagation (dt >= 0)
    static SolveStatus propagate(double mu, const State2& initial, double dt, State2& out_state);

    // Orbital invariants (return raw computed values; caller validates inputs)
    static double specific_energy(double mu, const State2& state);
    static double specific_angular_momentum_z(const State2& state);
    static double eccentricity(double mu, const State2& state);
    static double semi_major_axis(double mu, const State2& state);
    static double semi_latus_rectum(double mu, const State2& state);
    static double periapsis_radius(double mu, const State2& state);
    static double apoapsis_radius(double mu, const State2& state); // infinity for open orbits
};

// Internal helpers exposed for testing. Not part of the public API contract.
namespace detail {

// Angle normalization to (-pi, pi]
double wrap_angle(double a);

// Clamped inverse trig (inputs clamped to [-1, 1], never returns NaN)
double safe_acos(double x);
double safe_asin(double x);

// --- Elliptic anomaly conversions ---
// All assume 0 <= e < 1

double true_to_eccentric_anomaly(double nu, double e);
double eccentric_to_true_anomaly(double E, double e);
double eccentric_to_mean_anomaly(double E, double e);
SolveStatus mean_to_eccentric_anomaly(double M, double e, double& out_E,
                                      int max_iterations = 50);

// --- Hyperbolic anomaly conversions ---
// All assume e > 1

double true_to_hyperbolic_anomaly(double nu, double e);
double hyperbolic_to_true_anomaly(double H, double e);
double hyperbolic_to_mean_anomaly(double H, double e);
SolveStatus mean_to_hyperbolic_anomaly(double Mh, double e, double& out_H,
                                       int max_iterations = 50);

} // namespace detail

} // namespace brahe
