#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/trajectory_builder.h"
#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double kDay = 86400.0;
constexpr BodyId kEarth = 0;
constexpr BodyId kMoon = 1;

BodySystem earth_moon_system() {
    BodySystemBuilder b;

    BodyDef earth{};
    earth.id = kEarth;
    earth.parent_id = InvalidBody;
    earth.mu = 398600.4418;
    earth.radius = 6378.137;
    earth.soi_radius = 925000.0;
    b.add_body(earth);

    BodyDef moon{};
    moon.id = kMoon;
    moon.parent_id = kEarth;
    moon.mu = 4902.800066;
    moon.radius = 1737.4;
    moon.orbit_radius = 384400.0;
    moon.soi_radius = moon.orbit_radius * std::pow(moon.mu / earth.mu, 0.4);
    moon.angular_rate = std::sqrt(earth.mu / (moon.orbit_radius * moon.orbit_radius *
                                              moon.orbit_radius));
    moon.phase_at_epoch = 130.0 * M_PI / 180.0;
    b.add_body(moon);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);
    return sys;
}

BodySystem root_only(double mu, double radius, double soi_radius) {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = kEarth;
    root.parent_id = InvalidBody;
    root.mu = mu;
    root.radius = radius;
    root.soi_radius = soi_radius;
    b.add_body(root);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);
    return sys;
}

Tolerances preview_tolerances() {
    Tolerances tol{};
    tol.root_epsilon = 1e-3;
    tol.max_event_refine_iterations = 96;
    return tol;
}

}  // namespace

TEST_CASE("PreviewChain_FreeReturnImpactDoesNotDependOnEndTime",
          "[phase4][trajectory][central_radius]") {
    BodySystem sys = earth_moon_system();
    TrajectoryBuilder builder(sys, preview_tolerances());

    constexpr double expected_impact_day = 6.74580583190028;
    for (double end_day : {6.75, 6.8, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0, 14.0}) {
        PreviewRequest req;
        req.central_body = kEarth;
        req.start_time = 0.0;
        req.initial_state = State2{{6678.0, 0.0}, {0.0, 10.85}};
        req.end_time = end_day * kDay;
        req.max_segments = 8;

        Trajectory out;
        REQUIRE(builder.build_preview(req, out) == SolveStatus::Ok);
        REQUIRE(!out.segments.empty());
        REQUIRE(out.segments.back().end_reason == EventType::Impact);
        REQUIRE_THAT(out.segments.back().end_time / kDay,
                     WithinAbs(expected_impact_day, 1e-5));
    }
}

TEST_CASE("EventDetector_CentralImpactCrossingIsStableAcrossHorizons",
          "[phase4][event][central_radius]") {
    BodySystem sys = root_only(398600.4418, 6378.137, 925000.0);
    EventDetector detector(sys, preview_tolerances());

    constexpr double segment_start = 338564.2405093991;
    constexpr double expected_impact = 582837.6226361587;
    for (double end_day : {6.75, 6.8, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0, 14.0}) {
        EventSearchRequest req;
        req.central_body = kEarth;
        req.start_time = segment_start;
        req.initial_state = State2{{-353878.49528827524, 49475.43942851123},
                                   {0.598126549075723, -0.1798073770051678}};
        req.time_limit = end_day * kDay;

        PredictedEvent event;
        REQUIRE(detector.find_next_event(req, event) == SolveStatus::Ok);
        REQUIRE(event.type == EventType::Impact);
        REQUIRE_THAT(event.time, WithinAbs(expected_impact, 1.0));
    }
}

TEST_CASE("EventDetector_CentralSoiExitCrossingUsesApoapsisCandidate",
          "[phase4][event][central_radius]") {
    BodySystem sys = root_only(1.0, 0.1, 10.0);
    EventDetector detector(sys, preview_tolerances());

    // Ellipse with periapsis at r=5 and apoapsis just outside the SOI. The
    // central SOI exit is an apoapsis-side radius crossing and should not depend
    // on arbitrary coarse scan alignment.
    const double rp = 5.0;
    const double ra = 10.05;
    const double a = 0.5 * (rp + ra);
    const double vp = std::sqrt(1.0 * (2.0 / rp - 1.0 / a));

    for (double horizon : {40.0, 55.0, 70.0, 85.0, 100.0}) {
        EventSearchRequest req;
        req.central_body = kEarth;
        req.start_time = 0.0;
        req.initial_state = State2{{rp, 0.0}, {0.0, vp}};
        req.time_limit = horizon;

        PredictedEvent event;
        REQUIRE(detector.find_next_event(req, event) == SolveStatus::Ok);
        REQUIRE(event.type == EventType::SoiExit);
        REQUIRE_THAT(length(event.state.r), WithinAbs(10.0, 1e-3));
    }
}
