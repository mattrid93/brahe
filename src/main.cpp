#include <cstdio>

#include "brahe/body_system.h"
#include "brahe/vec2.h"

int main() {
    using namespace brahe;

    BodySystemBuilder builder;

    BodyDef star;
    star.id = 0;
    star.parent_id = InvalidBody;
    star.mu = 1.327e20;
    star.radius = 6.957e8;
    star.soi_radius = 1e12;
    builder.add_body(star);

    BodyDef planet;
    planet.id = 1;
    planet.parent_id = 0;
    planet.mu = 3.986e14;
    planet.radius = 6.371e6;
    planet.soi_radius = 9.246e8;
    planet.orbit_radius = 1.496e11;
    planet.angular_rate = 1.991e-7;
    planet.phase_at_epoch = 0.0;
    builder.add_body(planet);

    BodyDef moon;
    moon.id = 2;
    moon.parent_id = 1;
    moon.mu = 4.905e12;
    moon.radius = 1.737e6;
    moon.soi_radius = 6.616e7;
    moon.orbit_radius = 3.844e8;
    moon.angular_rate = 2.662e-6;
    moon.phase_at_epoch = 0.0;
    builder.add_body(moon);

    BodySystem sys;
    if (builder.build(sys) != SolveStatus::Ok) {
        std::printf("Failed to build body system\n");
        return 1;
    }

    std::printf("Body system built: %zu bodies, root=%u\n", sys.body_count(), sys.root_id());

    double t = 86400.0; // 1 day
    State2 planet_state = sys.state_in_root_frame(1, t);
    std::printf("Planet at t=%.0fs: pos=(%.3e, %.3e) vel=(%.3e, %.3e)\n", t, planet_state.r.x,
                planet_state.r.y, planet_state.v.x, planet_state.v.y);

    State2 moon_state = sys.state_in_root_frame(2, t);
    std::printf("Moon   at t=%.0fs: pos=(%.3e, %.3e) vel=(%.3e, %.3e)\n", t, moon_state.r.x,
                moon_state.r.y, moon_state.v.x, moon_state.v.y);

    return 0;
}
