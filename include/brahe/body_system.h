#pragma once

#include "brahe/types.h"

#include <cstddef>
#include <vector>

namespace brahe {

class BodySystem;

class BodySystemBuilder {
public:
    void add_body(const BodyDef& def);
    SolveStatus build(BodySystem& out) const;

private:
    std::vector<BodyDef> defs_;
};

class BodySystem {
public:
    const BodyDef* get_body(BodyId id) const;

    Vec2 position_in_parent(BodyId id, double t) const;
    Vec2 velocity_in_parent(BodyId id, double t) const;

    SolveStatus state_in_ancestor_frame(BodyId body, BodyId ancestor, double t,
                                        State2& out) const;
    State2 state_in_root_frame(BodyId body, double t) const;

    BodyId root_id() const;
    size_t body_count() const;
    size_t depth(BodyId id) const;

    const BodyId* children_begin(BodyId id) const;
    const BodyId* children_end(BodyId id) const;

private:
    friend class BodySystemBuilder;

    struct BodyEntry {
        BodyDef def;
        size_t parent_index;
        size_t depth;
        size_t children_start;
        size_t children_count;
    };

    std::vector<BodyEntry> bodies_;
    std::vector<BodyId> children_;
    size_t root_index_ = 0;

    size_t find_index(BodyId id) const;
    bool is_ancestor(size_t body_idx, size_t ancestor_idx) const;
};

} // namespace brahe
