#pragma once

#include <cstdint>
#include <limits>

namespace brahe {

using BodyId = uint32_t;
constexpr BodyId InvalidBody = ~BodyId{0};

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

struct State2 {
    Vec2 r;
    Vec2 v;
};

enum class ConicType { Ellipse, ParabolaLike, Hyperbola };

enum class EventType { None, Impact, SoiEntry, SoiExit, Burn, TimeLimit };

enum class SolveStatus {
    Ok,
    InvalidInput,
    NoConvergence,
    NoSolution,
    NumericalFailure,
    CapacityExceeded
};

enum class LambertDirection { CCW, CW };

enum class BurnFrame { Inertial, ProgradeRadial, ProgradeNormal };

struct Tolerances {
    double position_epsilon = 1e-3;
    double velocity_epsilon = 1e-6;
    double angle_epsilon = 1e-9;
    double time_epsilon = 1e-6;
    double root_epsilon = 1e-6;
    double lambert_residual_epsilon = 1e-10;
    double parabolic_eccentricity_band = 1e-6;
    int max_event_refine_iterations = 32;
};

struct BodyDef {
    BodyId id = InvalidBody;
    BodyId parent_id = InvalidBody;
    double mu = 0.0;
    double radius = 0.0;
    double soi_radius = 0.0;
    double orbit_radius = 0.0;
    double angular_rate = 0.0;
    double phase_at_epoch = 0.0;
};

} // namespace brahe
