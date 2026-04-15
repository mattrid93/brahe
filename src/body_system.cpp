#include "brahe/body_system.h"

#include "brahe/vec2.h"

#include <algorithm>
#include <cmath>
#include <queue>

namespace brahe {

// --- BodySystemBuilder ---

void BodySystemBuilder::add_body(const BodyDef& def) { defs_.push_back(def); }

SolveStatus BodySystemBuilder::build(BodySystem& out) const {
    if (defs_.empty()) {
        return SolveStatus::InvalidInput;
    }

    // Copy and sort by BodyId for canonical ordering
    std::vector<BodyDef> sorted = defs_;
    std::sort(sorted.begin(), sorted.end(),
              [](const BodyDef& a, const BodyDef& b) { return a.id < b.id; });

    // Check for duplicate IDs
    for (size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i].id == sorted[i - 1].id) {
            return SolveStatus::InvalidInput;
        }
    }

    size_t n = sorted.size();

    // Build id -> index map via binary search helper
    auto find_idx = [&](BodyId id) -> size_t {
        auto it = std::lower_bound(sorted.begin(), sorted.end(), id,
                                   [](const BodyDef& d, BodyId target) { return d.id < target; });
        if (it != sorted.end() && it->id == id) {
            return static_cast<size_t>(it - sorted.begin());
        }
        return SIZE_MAX;
    };

    // Find root and validate parent links
    size_t root_idx = SIZE_MAX;
    std::vector<size_t> parent_indices(n, SIZE_MAX);

    for (size_t i = 0; i < n; ++i) {
        if (sorted[i].parent_id == InvalidBody) {
            if (root_idx != SIZE_MAX) {
                return SolveStatus::InvalidInput; // multiple roots
            }
            root_idx = i;
        } else {
            size_t pi = find_idx(sorted[i].parent_id);
            if (pi == SIZE_MAX) {
                return SolveStatus::InvalidInput; // parent not found
            }
            parent_indices[i] = pi;
        }
    }

    if (root_idx == SIZE_MAX) {
        return SolveStatus::InvalidInput; // no root
    }

    // Validate constraints for non-root bodies
    for (size_t i = 0; i < n; ++i) {
        const auto& def = sorted[i];
        if (i == root_idx) continue;

        if (def.orbit_radius <= 0.0) {
            return SolveStatus::InvalidInput;
        }
        if (def.soi_radius <= def.radius) {
            return SolveStatus::InvalidInput;
        }
        if (def.soi_radius >= def.orbit_radius) {
            return SolveStatus::InvalidInput;
        }
    }

    // Compute depths via BFS from root
    std::vector<size_t> depths(n, SIZE_MAX);
    depths[root_idx] = 0;

    std::queue<size_t> bfs;
    bfs.push(root_idx);

    while (!bfs.empty()) {
        size_t cur = bfs.front();
        bfs.pop();

        for (size_t i = 0; i < n; ++i) {
            if (parent_indices[i] == cur && depths[i] == SIZE_MAX) {
                depths[i] = depths[cur] + 1;
                bfs.push(i);
            }
        }
    }

    // Any unvisited body means a disconnected cycle
    for (size_t i = 0; i < n; ++i) {
        if (depths[i] == SIZE_MAX) {
            return SolveStatus::InvalidInput;
        }
    }

    // Build child ranges
    // Count children per parent
    std::vector<size_t> child_counts(n, 0);
    for (size_t i = 0; i < n; ++i) {
        if (i != root_idx) {
            child_counts[parent_indices[i]]++;
        }
    }

    // Compute start offsets
    std::vector<size_t> child_starts(n);
    size_t total_children = 0;
    for (size_t i = 0; i < n; ++i) {
        child_starts[i] = total_children;
        total_children += child_counts[i];
    }

    // Fill child IDs (already in ascending BodyId order since sorted array is iterated in order)
    std::vector<BodyId> children(total_children);
    std::vector<size_t> fill_pos = child_starts;
    for (size_t i = 0; i < n; ++i) {
        if (i != root_idx) {
            size_t pi = parent_indices[i];
            children[fill_pos[pi]++] = sorted[i].id;
        }
    }

    // Build output
    out.bodies_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        out.bodies_[i].def = sorted[i];
        out.bodies_[i].parent_index = parent_indices[i];
        out.bodies_[i].depth = depths[i];
        out.bodies_[i].children_start = child_starts[i];
        out.bodies_[i].children_count = child_counts[i];
    }
    out.children_ = std::move(children);
    out.root_index_ = root_idx;

    return SolveStatus::Ok;
}

// --- BodySystem ---

size_t BodySystem::find_index(BodyId id) const {
    // Binary search on sorted bodies_
    size_t lo = 0;
    size_t hi = bodies_.size();
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (bodies_[mid].def.id < id) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo < bodies_.size() && bodies_[lo].def.id == id) {
        return lo;
    }
    return SIZE_MAX;
}

bool BodySystem::is_ancestor(size_t body_idx, size_t ancestor_idx) const {
    size_t cur = body_idx;
    while (cur != SIZE_MAX) {
        if (cur == ancestor_idx) return true;
        cur = bodies_[cur].parent_index;
    }
    return false;
}

const BodyDef* BodySystem::get_body(BodyId id) const {
    size_t idx = find_index(id);
    if (idx == SIZE_MAX) return nullptr;
    return &bodies_[idx].def;
}

Vec2 BodySystem::position_in_parent(BodyId id, double t) const {
    size_t idx = find_index(id);
    if (idx == SIZE_MAX || idx == root_index_) {
        return {0.0, 0.0};
    }
    const auto& def = bodies_[idx].def;
    double theta = def.phase_at_epoch + def.angular_rate * t;
    return {def.orbit_radius * std::cos(theta), def.orbit_radius * std::sin(theta)};
}

Vec2 BodySystem::velocity_in_parent(BodyId id, double t) const {
    size_t idx = find_index(id);
    if (idx == SIZE_MAX || idx == root_index_) {
        return {0.0, 0.0};
    }
    const auto& def = bodies_[idx].def;
    double theta = def.phase_at_epoch + def.angular_rate * t;
    double v_mag = def.orbit_radius * def.angular_rate;
    return {v_mag * (-std::sin(theta)), v_mag * std::cos(theta)};
}

SolveStatus BodySystem::state_in_ancestor_frame(BodyId body, BodyId ancestor, double t,
                                                 State2& out) const {
    size_t body_idx = find_index(body);
    if (body_idx == SIZE_MAX) return SolveStatus::InvalidInput;

    size_t ancestor_idx = find_index(ancestor);
    if (ancestor_idx == SIZE_MAX) return SolveStatus::InvalidInput;

    if (!is_ancestor(body_idx, ancestor_idx)) {
        return SolveStatus::InvalidInput;
    }

    // Accumulate position and velocity from body up to (but not including) ancestor
    Vec2 pos = {0.0, 0.0};
    Vec2 vel = {0.0, 0.0};
    size_t cur = body_idx;

    while (cur != ancestor_idx) {
        pos += position_in_parent(bodies_[cur].def.id, t);
        vel += velocity_in_parent(bodies_[cur].def.id, t);
        cur = bodies_[cur].parent_index;
    }

    out.r = pos;
    out.v = vel;
    return SolveStatus::Ok;
}

State2 BodySystem::state_in_root_frame(BodyId body, double t) const {
    State2 out{};
    state_in_ancestor_frame(body, bodies_[root_index_].def.id, t, out);
    return out;
}

BodyId BodySystem::root_id() const {
    if (bodies_.empty()) return InvalidBody;
    return bodies_[root_index_].def.id;
}

size_t BodySystem::body_count() const { return bodies_.size(); }

size_t BodySystem::depth(BodyId id) const {
    size_t idx = find_index(id);
    if (idx == SIZE_MAX) return SIZE_MAX;
    return bodies_[idx].depth;
}

const BodyId* BodySystem::children_begin(BodyId id) const {
    size_t idx = find_index(id);
    if (idx == SIZE_MAX || bodies_[idx].children_count == 0) return nullptr;
    return &children_[bodies_[idx].children_start];
}

const BodyId* BodySystem::children_end(BodyId id) const {
    size_t idx = find_index(id);
    if (idx == SIZE_MAX || bodies_[idx].children_count == 0) return nullptr;
    return &children_[bodies_[idx].children_start + bodies_[idx].children_count];
}

} // namespace brahe
