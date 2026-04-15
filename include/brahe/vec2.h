#pragma once

#include "brahe/types.h"

#include <cmath>

namespace brahe {

inline Vec2 operator+(const Vec2& a, const Vec2& b) { return {a.x + b.x, a.y + b.y}; }
inline Vec2 operator-(const Vec2& a, const Vec2& b) { return {a.x - b.x, a.y - b.y}; }
inline Vec2 operator-(const Vec2& v) { return {-v.x, -v.y}; }
inline Vec2 operator*(double s, const Vec2& v) { return {s * v.x, s * v.y}; }
inline Vec2 operator*(const Vec2& v, double s) { return {s * v.x, s * v.y}; }

inline Vec2& operator+=(Vec2& a, const Vec2& b) {
    a.x += b.x;
    a.y += b.y;
    return a;
}

inline double dot(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }
inline double length_squared(const Vec2& v) { return dot(v, v); }
inline double length(const Vec2& v) { return std::sqrt(length_squared(v)); }

inline bool is_finite(const Vec2& v) { return std::isfinite(v.x) && std::isfinite(v.y); }

inline bool normalize(const Vec2& v, Vec2& out) {
    double len = length(v);
    if (len < 1e-15) {
        out = {0.0, 0.0};
        return false;
    }
    out = {v.x / len, v.y / len};
    return true;
}

} // namespace brahe
