#pragma once

#include "brahe/types.h"

namespace brahe {

struct ConicElements2D {
    ConicType type = ConicType::Ellipse;

    double mu = 0.0;
    double semi_major_axis = 0.0;
    double eccentricity = 0.0;
    double semi_latus_rectum = 0.0;
    double angular_momentum_z = 0.0; // h = x*vy - y*vx; sign encodes orbit direction
    double argument_of_periapsis = 0.0; // (-pi, pi]
    double true_anomaly = 0.0;         // (-pi, pi]

    double mean_anomaly = 0.0;      // defined for elliptic
    double eccentric_anomaly = 0.0;  // defined for elliptic
    double hyperbolic_anomaly = 0.0; // defined for hyperbolic
};

struct ConicState {
    BodyId central_body = InvalidBody;
    double epoch = 0.0;
    State2 state_at_epoch;
    ConicType type = ConicType::Ellipse;
};

} // namespace brahe
