#include <catch2/catch_test_macros.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/types.h"

#include <cmath>

using namespace brahe;

namespace {

constexpr BodyId kPlanet = 0;
constexpr BodyId kMoon = 1;

BodySystem single_moon_system() {
    BodyDef planet;
    planet.id = kPlanet;
    planet.parent_id = InvalidBody;
    planet.mu = 126686534.0;
    planet.radius = 71492.0;
    planet.soi_radius = 48000000.0;

    BodyDef moon;
    moon.id = kMoon;
    moon.parent_id = kPlanet;
    moon.mu = 3200.0;
    moon.radius = 1400.0;
    moon.soi_radius = 24000.0;
    moon.orbit_radius = 310000.0;
    moon.angular_rate = std::sqrt(planet.mu / (moon.orbit_radius * moon.orbit_radius *
                                               moon.orbit_radius));
    moon.phase_at_epoch = 0.20;

    BodySystemBuilder builder;
    builder.add_body(planet);
    builder.add_body(moon);

    BodySystem sys;
    REQUIRE(builder.build(sys) == SolveStatus::Ok);
    return sys;
}

}  // namespace

TEST_CASE("EventDetector_SkipsCoarseScanWhenChildSoiRadiallyUnreachable",
          "[event][performance]") {
    BodySystem sys = single_moon_system();
    Tolerances tol{};
    tol.root_epsilon = 1e-3;
    tol.max_event_refine_iterations = 96;
    EventDetectorDiagnostics diagnostics{};
    EventDetector detector(sys, tol, &diagnostics);

    EventSearchRequest req;
    req.central_body = kPlanet;
    req.start_time = 0.0;
    req.time_limit = 300.0 * 86400.0;
    req.initial_state = State2{{242598.819617 * std::cos(1.352099192),
                                242598.819617 * std::sin(1.352099192)},
                               {-21.177021862 * std::sin(1.352099192),
                                21.177021862 * std::cos(1.352099192)}};

    PredictedEvent event;
    REQUIRE(detector.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::TimeLimit);
    REQUIRE(diagnostics.scan_steps == 0);
}
