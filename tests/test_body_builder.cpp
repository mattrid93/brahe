#include <catch2/catch_test_macros.hpp>

#include "brahe/body_system.h"

#include <cmath>
#include <limits>

using namespace brahe;

// Helper: minimal valid star-planet-moon system
static void add_star_planet_moon(BodySystemBuilder& b) {
    BodyDef star;
    star.id = 0;
    star.parent_id = InvalidBody;
    star.mu = 1.0e10;
    star.radius = 100.0;
    star.soi_radius = 1e9;

    BodyDef planet;
    planet.id = 1;
    planet.parent_id = 0;
    planet.mu = 1.0e6;
    planet.radius = 10.0;
    planet.soi_radius = 1000.0;
    planet.orbit_radius = 10000.0;
    planet.angular_rate = 0.001;
    planet.phase_at_epoch = 0.0;

    BodyDef moon;
    moon.id = 2;
    moon.parent_id = 1;
    moon.mu = 1.0e3;
    moon.radius = 1.0;
    moon.soi_radius = 50.0;
    moon.orbit_radius = 500.0;
    moon.angular_rate = 0.01;
    moon.phase_at_epoch = 0.0;

    b.add_body(star);
    b.add_body(planet);
    b.add_body(moon);
}

static BodyDef valid_root() {
    BodyDef root;
    root.id = 0;
    root.parent_id = InvalidBody;
    root.mu = 1.0e10;
    root.radius = 100.0;
    root.soi_radius = 1e9;
    return root;
}

static BodyDef valid_child() {
    BodyDef child;
    child.id = 1;
    child.parent_id = 0;
    child.mu = 1.0e6;
    child.radius = 10.0;
    child.soi_radius = 1000.0;
    child.orbit_radius = 10000.0;
    child.angular_rate = 0.001;
    child.phase_at_epoch = 0.0;
    return child;
}

// --- B. Builder and hierarchy validation tests ---

TEST_CASE("Build_Fails_OnDuplicateBodyId", "[builder]") {
    BodySystemBuilder b;
    BodyDef d1;
    d1.id = 0;
    d1.parent_id = InvalidBody;
    d1.mu = 1.0;
    d1.radius = 1.0;
    d1.soi_radius = 10.0;

    BodyDef d2 = d1; // same id
    b.add_body(d1);
    b.add_body(d2);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
}

TEST_CASE("Build_Fails_WhenParentMissing", "[builder]") {
    BodySystemBuilder b;
    BodyDef root;
    root.id = 0;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 1.0;
    root.soi_radius = 10.0;

    BodyDef child;
    child.id = 1;
    child.parent_id = 99; // nonexistent
    child.mu = 1.0;
    child.radius = 1.0;
    child.soi_radius = 5.0;
    child.orbit_radius = 100.0;
    child.angular_rate = 0.01;

    b.add_body(root);
    b.add_body(child);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
}

TEST_CASE("Build_Fails_WhenNoRootExists", "[builder]") {
    BodySystemBuilder b;
    BodyDef d;
    d.id = 0;
    d.parent_id = 1;
    d.mu = 1.0;
    d.radius = 1.0;
    d.soi_radius = 5.0;
    d.orbit_radius = 100.0;
    d.angular_rate = 0.01;

    b.add_body(d);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
}

TEST_CASE("Build_Fails_WhenMultipleRootsExist", "[builder]") {
    BodySystemBuilder b;
    BodyDef r1;
    r1.id = 0;
    r1.parent_id = InvalidBody;
    r1.mu = 1.0;
    r1.radius = 1.0;
    r1.soi_radius = 10.0;

    BodyDef r2;
    r2.id = 1;
    r2.parent_id = InvalidBody;
    r2.mu = 1.0;
    r2.radius = 1.0;
    r2.soi_radius = 10.0;

    b.add_body(r1);
    b.add_body(r2);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
}

TEST_CASE("Build_Fails_OnCycle", "[builder]") {
    BodySystemBuilder b;

    // Root exists, but A and B form a cycle disconnected from root
    BodyDef root;
    root.id = 0;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 1.0;
    root.soi_radius = 10.0;

    BodyDef a;
    a.id = 1;
    a.parent_id = 2;
    a.mu = 1.0;
    a.radius = 1.0;
    a.soi_radius = 5.0;
    a.orbit_radius = 100.0;
    a.angular_rate = 0.01;

    BodyDef bb;
    bb.id = 2;
    bb.parent_id = 1;
    bb.mu = 1.0;
    bb.radius = 1.0;
    bb.soi_radius = 5.0;
    bb.orbit_radius = 100.0;
    bb.angular_rate = 0.01;

    b.add_body(root);
    b.add_body(a);
    b.add_body(bb);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
}

TEST_CASE("Build_Fails_WhenNonRootOrbitRadiusNotPositive", "[builder]") {
    BodySystemBuilder b;
    BodyDef root;
    root.id = 0;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 1.0;
    root.soi_radius = 10.0;

    BodyDef child;
    child.id = 1;
    child.parent_id = 0;
    child.mu = 1.0;
    child.radius = 1.0;
    child.soi_radius = 5.0;
    child.orbit_radius = 0.0; // invalid
    child.angular_rate = 0.01;

    b.add_body(root);
    b.add_body(child);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
}

TEST_CASE("Build_Fails_WhenSoiNotLargerThanRadius", "[builder]") {
    BodySystemBuilder b;
    BodyDef root;
    root.id = 0;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 1.0;
    root.soi_radius = 10.0;

    BodyDef child;
    child.id = 1;
    child.parent_id = 0;
    child.mu = 1.0;
    child.radius = 10.0;
    child.soi_radius = 10.0; // not larger than radius
    child.orbit_radius = 100.0;
    child.angular_rate = 0.01;

    b.add_body(root);
    b.add_body(child);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
}

TEST_CASE("Build_Fails_WhenChildSoiDoesNotFitOrbit", "[builder]") {
    BodySystemBuilder b;
    BodyDef root;
    root.id = 0;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 1.0;
    root.soi_radius = 10.0;

    BodyDef child;
    child.id = 1;
    child.parent_id = 0;
    child.mu = 1.0;
    child.radius = 1.0;
    child.soi_radius = 100.0; // >= orbit_radius
    child.orbit_radius = 100.0;
    child.angular_rate = 0.01;

    b.add_body(root);
    b.add_body(child);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
}

TEST_CASE("Build_Fails_WhenRootSoiNotLargerThanRadius", "[builder]") {
    BodySystemBuilder b;
    BodyDef root = valid_root();
    root.soi_radius = root.radius;

    b.add_body(root);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
}

TEST_CASE("Build_Fails_WhenBodyNumericFieldsAreInvalid", "[builder]") {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    auto require_invalid = [](BodyDef root, BodyDef child) {
        BodySystemBuilder b;
        b.add_body(root);
        b.add_body(child);
        BodySystem sys;
        REQUIRE(b.build(sys) == SolveStatus::InvalidInput);
    };

    SECTION("root mu must be positive and finite") {
        BodyDef root = valid_root();
        root.mu = 0.0;
        require_invalid(root, valid_child());

        root = valid_root();
        root.mu = nan;
        require_invalid(root, valid_child());
    }

    SECTION("root radius and SOI must be finite") {
        BodyDef root = valid_root();
        root.radius = -1.0;
        require_invalid(root, valid_child());

        root = valid_root();
        root.soi_radius = inf;
        require_invalid(root, valid_child());
    }

    SECTION("child core fields must be positive and finite") {
        BodyDef child = valid_child();
        child.mu = -1.0;
        require_invalid(valid_root(), child);

        child = valid_child();
        child.radius = nan;
        require_invalid(valid_root(), child);

        child = valid_child();
        child.soi_radius = inf;
        require_invalid(valid_root(), child);
    }

    SECTION("child orbit fields must be finite") {
        BodyDef child = valid_child();
        child.orbit_radius = nan;
        require_invalid(valid_root(), child);

        child = valid_child();
        child.angular_rate = inf;
        require_invalid(valid_root(), child);

        child = valid_child();
        child.phase_at_epoch = nan;
        require_invalid(valid_root(), child);
    }
}

TEST_CASE("Build_Fails_WhenSiblingChildSpheresMayOverlap", "[builder]") {
    BodyDef root = valid_root();

    auto require_invalid = [&](BodyDef a, BodyDef b) {
        BodySystemBuilder builder;
        builder.add_body(root);
        builder.add_body(a);
        builder.add_body(b);
        BodySystem sys;
        REQUIRE(builder.build(sys) == SolveStatus::InvalidInput);
    };

    SECTION("different orbit radii can bring SOIs into contact") {
        BodyDef a = valid_child();
        a.id = 1;
        a.orbit_radius = 1000.0;
        a.soi_radius = 60.0;
        a.radius = 10.0;

        BodyDef b = valid_child();
        b.id = 2;
        b.orbit_radius = 1090.0;
        b.soi_radius = 40.0;
        b.radius = 10.0;

        require_invalid(a, b);
    }

    SECTION("same angular rate and phase can make physical radii overlap") {
        BodyDef a = valid_child();
        a.id = 1;
        a.orbit_radius = 1000.0;
        a.angular_rate = 0.01;
        a.phase_at_epoch = 0.0;
        a.soi_radius = 20.0;
        a.radius = 10.0;

        BodyDef b = valid_child();
        b.id = 2;
        b.orbit_radius = 1000.0;
        b.angular_rate = 0.01;
        b.phase_at_epoch = 0.0;
        b.soi_radius = 20.0;
        b.radius = 10.0;

        require_invalid(a, b);
    }
}

TEST_CASE("Build_Succeeds_OnValidStarPlanetMoonTree", "[builder]") {
    BodySystemBuilder b;
    add_star_planet_moon(b);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);
    REQUIRE(sys.body_count() == 3);
    REQUIRE(sys.root_id() == 0);
}

TEST_CASE("GetBody_ReturnsNullForUnknownId", "[builder]") {
    BodySystemBuilder b;
    add_star_planet_moon(b);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);
    REQUIRE(sys.get_body(999) == nullptr);
}

TEST_CASE("GetBody_ReturnsStableDefinitionForKnownId", "[builder]") {
    BodySystemBuilder b;
    add_star_planet_moon(b);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);

    const BodyDef* planet = sys.get_body(1);
    REQUIRE(planet != nullptr);
    REQUIRE(planet->id == 1);
    REQUIRE(planet->parent_id == 0);
    REQUIRE(planet->mu == 1.0e6);
    REQUIRE(planet->orbit_radius == 10000.0);
    REQUIRE(planet->angular_rate == 0.001);
}

// --- C. Deterministic ordering and cached topology tests ---

TEST_CASE("Build_CanonicalizesStorageIndependentOfInsertionOrder", "[builder]") {
    // Build system in two different orders
    auto build_system = [](std::vector<int> order) {
        BodyDef star;
        star.id = 0;
        star.parent_id = InvalidBody;
        star.mu = 1.0e10;
        star.radius = 100.0;
        star.soi_radius = 1e9;

        BodyDef p1;
        p1.id = 1;
        p1.parent_id = 0;
        p1.mu = 1.0e6;
        p1.radius = 10.0;
        p1.soi_radius = 1000.0;
        p1.orbit_radius = 10000.0;
        p1.angular_rate = 0.001;
        p1.phase_at_epoch = 0.5;

        BodyDef p2;
        p2.id = 2;
        p2.parent_id = 0;
        p2.mu = 2.0e6;
        p2.radius = 20.0;
        p2.soi_radius = 2000.0;
        p2.orbit_radius = 20000.0;
        p2.angular_rate = 0.0005;
        p2.phase_at_epoch = 1.0;

        BodyDef defs[] = {star, p1, p2};
        BodySystemBuilder b;
        for (int i : order) b.add_body(defs[i]);

        BodySystem sys;
        b.build(sys);
        return sys;
    };

    BodySystem a = build_system({0, 1, 2});
    BodySystem b = build_system({2, 0, 1});

    // Same queries should produce identical results
    double t = 5.0;
    for (BodyId id : {BodyId(0), BodyId(1), BodyId(2)}) {
        Vec2 pa = a.position_in_parent(id, t);
        Vec2 pb = b.position_in_parent(id, t);
        REQUIRE(pa.x == pb.x);
        REQUIRE(pa.y == pb.y);

        Vec2 va = a.velocity_in_parent(id, t);
        Vec2 vb = b.velocity_in_parent(id, t);
        REQUIRE(va.x == vb.x);
        REQUIRE(va.y == vb.y);
    }

    State2 sa = a.state_in_root_frame(2, t);
    State2 sb = b.state_in_root_frame(2, t);
    REQUIRE(sa.r.x == sb.r.x);
    REQUIRE(sa.r.y == sb.r.y);
    REQUIRE(sa.v.x == sb.v.x);
    REQUIRE(sa.v.y == sb.v.y);
}

TEST_CASE("SiblingTraversalOrder_IsAscendingBodyId", "[builder]") {
    BodySystemBuilder b;

    BodyDef star;
    star.id = 0;
    star.parent_id = InvalidBody;
    star.mu = 1.0e10;
    star.radius = 100.0;
    star.soi_radius = 1e9;

    // Add children in reverse order
    BodyDef p3;
    p3.id = 30;
    p3.parent_id = 0;
    p3.mu = 1.0;
    p3.radius = 1.0;
    p3.soi_radius = 5.0;
    p3.orbit_radius = 300.0;
    p3.angular_rate = 0.01;

    BodyDef p1;
    p1.id = 10;
    p1.parent_id = 0;
    p1.mu = 1.0;
    p1.radius = 1.0;
    p1.soi_radius = 5.0;
    p1.orbit_radius = 100.0;
    p1.angular_rate = 0.01;

    BodyDef p2;
    p2.id = 20;
    p2.parent_id = 0;
    p2.mu = 1.0;
    p2.radius = 1.0;
    p2.soi_radius = 5.0;
    p2.orbit_radius = 200.0;
    p2.angular_rate = 0.01;

    b.add_body(p3);
    b.add_body(star);
    b.add_body(p1);
    b.add_body(p2);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);

    const BodyId* begin = sys.children_begin(0);
    const BodyId* end = sys.children_end(0);
    REQUIRE(begin != nullptr);
    REQUIRE(end != nullptr);

    std::vector<BodyId> children(begin, end);
    REQUIRE(children.size() == 3);
    REQUIRE(children[0] == 10);
    REQUIRE(children[1] == 20);
    REQUIRE(children[2] == 30);
}

TEST_CASE("ChildTraversalRange_ForLeafBodyIsStableEmptyRange", "[builder]") {
    BodySystemBuilder b;
    add_star_planet_moon(b);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);

    const BodyId* begin = sys.children_begin(2);
    const BodyId* end = sys.children_end(2);

    REQUIRE(begin != nullptr);
    REQUIRE(end != nullptr);
    REQUIRE(begin == end);
    REQUIRE(static_cast<size_t>(end - begin) == 0);
}

TEST_CASE("DepthCache_IsCorrectForNestedTree", "[builder]") {
    BodySystemBuilder b;
    add_star_planet_moon(b);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);

    REQUIRE(sys.depth(0) == 0); // root
    REQUIRE(sys.depth(1) == 1); // planet
    REQUIRE(sys.depth(2) == 2); // moon
}

TEST_CASE("AncestorQuery_RejectsNonAncestor", "[builder]") {
    BodySystemBuilder b;

    BodyDef star;
    star.id = 0;
    star.parent_id = InvalidBody;
    star.mu = 1.0e10;
    star.radius = 100.0;
    star.soi_radius = 1e9;

    BodyDef p1;
    p1.id = 1;
    p1.parent_id = 0;
    p1.mu = 1.0;
    p1.radius = 1.0;
    p1.soi_radius = 5.0;
    p1.orbit_radius = 100.0;
    p1.angular_rate = 0.01;

    BodyDef p2;
    p2.id = 2;
    p2.parent_id = 0;
    p2.mu = 1.0;
    p2.radius = 1.0;
    p2.soi_radius = 5.0;
    p2.orbit_radius = 200.0;
    p2.angular_rate = 0.01;

    b.add_body(star);
    b.add_body(p1);
    b.add_body(p2);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);

    // p2 is not an ancestor of p1
    State2 out;
    REQUIRE(sys.state_in_ancestor_frame(1, 2, 0.0, out) == SolveStatus::InvalidInput);
}
