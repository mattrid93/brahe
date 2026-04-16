#include <catch2/catch_test_macros.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"

#include <cmath>
#include <cstdint>
#include <random>

using namespace brahe;

// --- Section 7.12: Fuzz tests ---
// Deterministic PRNG so failures are reproducible.

namespace {

BodySystem build_root_only(double mu, double radius, double soi) {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = mu;
    root.radius = radius;
    root.soi_radius = soi;
    b.add_body(root);
    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

BodySystem build_with_children(int num_children, double child_soi, double orbit_radius,
                               std::mt19937* rng = nullptr) {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 0.5;
    root.soi_radius = 1000.0;
    b.add_body(root);

    // When rng is supplied, draw per-child angular_rate from [-0.2, 0.2] to exercise the
    // time-varying child ephemeris path (both prograde and retrograde).
    std::uniform_real_distribution<double> rate_dist(-0.2, 0.2);

    for (int i = 0; i < num_children; ++i) {
        BodyDef ch{};
        ch.id = static_cast<BodyId>(10 + i);
        ch.parent_id = 1;
        ch.mu = 0.001;
        ch.radius = 0.1;
        ch.soi_radius = child_soi;
        ch.orbit_radius = orbit_radius;
        ch.angular_rate = rng ? rate_dist(*rng) : 0.0;
        ch.phase_at_epoch = (2.0 * M_PI * i) / num_children;
        b.add_body(ch);
    }
    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

bool is_valid_output(const PredictedEvent& ev) {
    if (!std::isfinite(ev.time)) return false;
    if (!std::isfinite(ev.state.r.x) || !std::isfinite(ev.state.r.y)) return false;
    if (!std::isfinite(ev.state.v.x) || !std::isfinite(ev.state.v.y)) return false;
    return true;
}

}  // namespace

TEST_CASE("RandomValidSystemsAndStates_NoCrashesNoNaNs", "[phase3][fuzz]") {
    BodySystem sys = build_root_only(1.0, 0.5, 100.0);
    EventDetector det(sys);

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pos_dist(-20.0, 20.0);
    std::uniform_real_distribution<double> vel_dist(-2.0, 2.0);
    std::uniform_real_distribution<double> horizon_dist(1.0, 50.0);

    for (int i = 0; i < 200; ++i) {
        EventSearchRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        req.initial_state =
            State2{{pos_dist(rng), pos_dist(rng)}, {vel_dist(rng), vel_dist(rng)}};
        req.time_limit = horizon_dist(rng);

        PredictedEvent ev;
        auto status = det.find_next_event(req, ev);
        // Status must be a known value; output must not contain NaN even on failure.
        REQUIRE((status == SolveStatus::Ok || status == SolveStatus::InvalidInput ||
                 status == SolveStatus::NoConvergence ||
                 status == SolveStatus::NumericalFailure ||
                 status == SolveStatus::NoSolution));
        REQUIRE(is_valid_output(ev));
    }
}

TEST_CASE("RandomHighEccentricityCases_NoSkippedEventsWithinLimit", "[phase3][fuzz]") {
    BodySystem sys = build_root_only(1.0, 0.5, 50.0);
    EventDetector det(sys);

    std::mt19937 rng(67890);
    std::uniform_real_distribution<double> rp_dist(0.6, 1.5);    // above body radius 0.5
    std::uniform_real_distribution<double> ra_dist(20.0, 40.0);  // well inside SOI

    for (int i = 0; i < 100; ++i) {
        double rp = rp_dist(rng);
        double ra = ra_dist(rng);
        double a = 0.5 * (rp + ra);
        double e = (ra - rp) / (ra + rp);
        double vp = std::sqrt(1.0 * (1.0 + e) / (a * (1.0 - e)));

        EventSearchRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        req.initial_state = State2{{rp, 0.0}, {0.0, vp}};  // at periapsis
        req.time_limit = 100.0;

        PredictedEvent ev;
        auto status = det.find_next_event(req, ev);
        REQUIRE(is_valid_output(ev));
        // All these orbits are bound with periapsis above body radius and apoapsis
        // below SOI radius -- no actual event should occur -> TimeLimit.
        if (status == SolveStatus::Ok) {
            REQUIRE(ev.type == EventType::TimeLimit);
        }
    }
}

TEST_CASE("RandomSiblingChildSystems_EventOrderingAlwaysReturnsValidBodyIds", "[phase3][fuzz]") {
    std::mt19937 rng(2468);
    std::uniform_real_distribution<double> vel_dist(0.5, 5.0);
    std::uniform_int_distribution<int> nchild_dist(2, 4);

    for (int i = 0; i < 50; ++i) {
        int n = nchild_dist(rng);
        // Pass rng so the fuzz draws moving-moon angular rates as well.
        BodySystem sys = build_with_children(n, 5.0, 50.0, &rng);
        EventDetector det(sys);

        EventSearchRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        double vx = vel_dist(rng);
        req.initial_state = State2{{30.0, 0.0}, {vx, 0.0}};
        req.time_limit = 20.0;

        PredictedEvent ev;
        auto status = det.find_next_event(req, ev);
        REQUIRE(is_valid_output(ev));
        if (status == SolveStatus::Ok && ev.type == EventType::SoiEntry) {
            // to_body must be one of the configured child IDs (10..10+n-1)
            REQUIRE(ev.to_body >= 10);
            REQUIRE(ev.to_body < static_cast<BodyId>(10 + n));
        }
    }
}

TEST_CASE("RandomBoundaryStarts_DoNotInfiniteLoop", "[phase3][fuzz]") {
    BodySystem sys = build_root_only(1.0, 1.0, 10.0);
    EventDetector det(sys);

    std::mt19937 rng(13579);
    std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> vel_dist(-1.5, 1.5);

    for (int i = 0; i < 100; ++i) {
        double theta = angle_dist(rng);
        // Start exactly on a boundary (impact or SOI).
        double boundary = (i % 2 == 0) ? 1.0 : 10.0;
        Vec2 r = {boundary * std::cos(theta), boundary * std::sin(theta)};
        Vec2 v = {vel_dist(rng), vel_dist(rng)};

        EventSearchRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        req.initial_state = State2{r, v};
        req.time_limit = 20.0;

        PredictedEvent ev;
        auto status = det.find_next_event(req, ev);
        REQUIRE(is_valid_output(ev));
        (void)status;  // may succeed or fail; key property is completion without hang/NaN
    }
}
