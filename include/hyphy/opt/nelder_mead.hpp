#pragma once

#include "hyphy/core/types.hpp"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <functional>

namespace hyphy::opt {

using namespace hyphy::core;

// 1D Brent's method / Golden section search for bounded minimization
class Brent1D {
public:
    static std::pair<Scalar, Scalar> minimize(
        std::function<Scalar(Scalar)> f,
        Scalar ax, Scalar bx, Scalar cx,
        Scalar tol = 1e-4, int max_iter = 60
    ) {
        const Scalar CGOLD = 0.3819660;
        const Scalar ZEPS = 1e-10;

        Scalar a = std::min(ax, cx);
        Scalar b = std::max(ax, cx);
        Scalar v = bx, w = v, x = v;
        Scalar e = 0.0;
        Scalar fx = f(x), fv = fx, fw = fx;

        for (int iter = 0; iter < max_iter; ++iter) {
            Scalar xm = 0.5 * (a + b);
            Scalar tol1 = tol * std::abs(x) + ZEPS;
            Scalar tol2 = 2.0 * tol1;

            if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                return {x, fx};
            }

            Scalar d = 0.0;
            if (std::abs(e) > tol1) {
                Scalar r = (x - w) * (fx - fv);
                Scalar q = (x - v) * (fx - fw);
                Scalar p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if (q > 0.0) p = -p;
                q = std::abs(q);
                Scalar etemp = e;
                e = d;

                if (std::abs(p) >= std::abs(0.5 * q * etemp) ||
                    p <= q * (a - x) || p >= q * (b - x)) {
                    e = (x >= xm) ? a - x : b - x;
                    d = CGOLD * e;
                } else {
                    d = p / q;
                    Scalar u = x + d;
                    if (u - a < tol2 || b - u < tol2) {
                        d = (xm - x >= 0.0) ? std::abs(tol1) : -std::abs(tol1);
                    }
                }
            } else {
                e = (x >= xm) ? a - x : b - x;
                d = CGOLD * e;
            }

            Scalar u = (std::abs(d) >= tol1) ? x + d : x + ((d >= 0.0) ? std::abs(tol1) : -std::abs(tol1));
            Scalar fu = f(u);

            if (fu <= fx) {
                if (u >= x) a = x; else b = x;
                v = w; fv = fw;
                w = x; fw = fx;
                x = u; fx = fu;
            } else {
                if (u < x) a = u; else b = u;
                if (fu <= fw || w == x) {
                    v = w; fv = fw;
                    w = u; fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    v = u; fv = fu;
                }
            }
        }

        return {x, fx};
    }
};

// 2D Nelder-Mead simplex optimizer with bound constraints
class NelderMead2D {
public:
    struct Point {
        Scalar a = 0.0;
        Scalar b = 0.0;
        Scalar f = 0.0;
    };

    static Point minimize(
        std::function<Scalar(Scalar, Scalar)> f,
        Scalar init_a, Scalar init_b,
        Scalar lb = 1e-4, Scalar ub = 1000.0,
        Scalar tol = 1e-4, int max_iter = 100
    ) {
        auto clamp = [&](Scalar v) {
            return std::clamp(v, lb, ub);
        };

        auto eval = [&](Scalar a, Scalar b) -> Point {
            Scalar ca = clamp(a);
            Scalar cb = clamp(b);
            return {ca, cb, f(ca, cb)};
        };

        // Construct initial simplex around (init_a, init_b)
        Scalar step_a = std::max(0.1, 0.2 * init_a);
        Scalar step_b = std::max(0.1, 0.2 * init_b);

        std::array<Point, 3> s = {
            eval(init_a, init_b),
            eval(init_a + step_a, init_b),
            eval(init_a, init_b + step_b)
        };

        const Scalar alpha = 1.0; // reflection
        const Scalar gamma = 2.0; // expansion
        const Scalar rho = 0.5;   // contraction
        const Scalar sigma = 0.5; // shrink

        for (int iter = 0; iter < max_iter; ++iter) {
            // Sort: s[0] best, s[2] worst
            std::sort(s.begin(), s.end(), [](const Point& p1, const Point& p2) {
                return p1.f < p2.f;
            });

            if (std::abs(s[2].f - s[0].f) < tol &&
                std::abs(s[2].a - s[0].a) < tol &&
                std::abs(s[2].b - s[0].b) < tol) {
                break;
            }

            // Centroid of s[0] and s[1]
            Scalar c_a = 0.5 * (s[0].a + s[1].a);
            Scalar c_b = 0.5 * (s[0].b + s[1].b);

            // Reflection
            Point r = eval(c_a + alpha * (c_a - s[2].a), c_b + alpha * (c_b - s[2].b));

            if (r.f < s[0].f) {
                // Expansion
                Point e = eval(c_a + gamma * (r.a - c_a), c_b + gamma * (r.b - c_b));
                s[2] = (e.f < r.f) ? e : r;
            } else if (r.f < s[1].f) {
                s[2] = r;
            } else {
                // Contraction
                bool outside = (r.f < s[2].f);
                Point c = outside
                    ? eval(c_a + rho * (r.a - c_a), c_b + rho * (r.b - c_b))
                    : eval(c_a + rho * (s[2].a - c_a), c_b + rho * (s[2].b - c_b));

                Scalar ref_f = outside ? r.f : s[2].f;
                if (c.f < ref_f) {
                    s[2] = c;
                } else {
                    // Shrink
                    s[1] = eval(s[0].a + sigma * (s[1].a - s[0].a), s[0].b + sigma * (s[1].b - s[0].b));
                    s[2] = eval(s[0].a + sigma * (s[2].a - s[0].a), s[0].b + sigma * (s[2].b - s[0].b));
                }
            }
        }

        std::sort(s.begin(), s.end(), [](const Point& p1, const Point& p2) {
            return p1.f < p2.f;
        });
        return s[0];
    }
};

// General N-dimensional Nelder-Mead simplex optimizer with box constraints and initial grid
template <size_t Dim>
class NelderMeadND {
public:
    struct Point {
        std::array<Scalar, Dim> x{};
        Scalar f = 0.0;
    };

    static Point minimize(
        std::function<Scalar(const std::array<Scalar, Dim>&)> func,
        const std::array<Scalar, Dim>& init_x,
        const std::array<Scalar, Dim>& lb,
        const std::array<Scalar, Dim>& ub,
        Scalar tol = 1e-4,
        int max_iter = 150,
        const std::vector<std::array<Scalar, Dim>>& additional_grid_points = {}
    ) {
        auto clamp = [&](const std::array<Scalar, Dim>& v) -> std::array<Scalar, Dim> {
            std::array<Scalar, Dim> cv;
            for (size_t d = 0; d < Dim; ++d) {
                cv[d] = std::clamp(v[d], lb[d], ub[d]);
            }
            return cv;
        };

        auto eval = [&](const std::array<Scalar, Dim>& v) -> Point {
            auto cv = clamp(v);
            return {cv, func(cv)};
        };

        // Determine best starting point among init_x and additional_grid_points
        Point best_start = eval(init_x);
        for (const auto& pt : additional_grid_points) {
            Point p = eval(pt);
            if (p.f < best_start.f) {
                best_start = p;
            }
        }

        // Build simplex with Dim + 1 points around best_start
        constexpr size_t N_verts = Dim + 1;
        std::array<Point, N_verts> s;
        s[0] = best_start;

        for (size_t d = 0; d < Dim; ++d) {
            std::array<Scalar, Dim> perturbed = best_start.x;
            Scalar step = std::max(0.05, 0.2 * std::abs(best_start.x[d]));
            if (best_start.x[d] + step <= ub[d]) {
                perturbed[d] += step;
            } else {
                perturbed[d] -= step;
            }
            s[d + 1] = eval(perturbed);
        }

        const Scalar alpha = 1.0; // reflection
        const Scalar gamma = 2.0; // expansion
        const Scalar rho   = 0.5; // contraction
        const Scalar sigma = 0.5; // shrink

        for (int iter = 0; iter < max_iter; ++iter) {
            // Sort vertices so s[0] is best, s[Dim] is worst
            std::sort(s.begin(), s.end(), [](const Point& p1, const Point& p2) {
                return p1.f < p2.f;
            });

            // Convergence check
            Scalar max_x_diff = 0.0;
            for (size_t d = 0; d < Dim; ++d) {
                max_x_diff = std::max(max_x_diff, std::abs(s[Dim].x[d] - s[0].x[d]));
            }
            if (std::abs(s[Dim].f - s[0].f) < tol && max_x_diff < tol) {
                break;
            }

            // Centroid of best Dim vertices
            std::array<Scalar, Dim> centroid{};
            for (size_t i = 0; i < Dim; ++i) {
                for (size_t d = 0; d < Dim; ++d) {
                    centroid[d] += s[i].x[d];
                }
            }
            for (size_t d = 0; d < Dim; ++d) {
                centroid[d] /= static_cast<Scalar>(Dim);
            }

            // Reflection: x_r = centroid + alpha * (centroid - s[Dim].x)
            std::array<Scalar, Dim> x_r;
            for (size_t d = 0; d < Dim; ++d) {
                x_r[d] = centroid[d] + alpha * (centroid[d] - s[Dim].x[d]);
            }
            Point r = eval(x_r);

            if (r.f < s[0].f) {
                // Expansion: x_e = centroid + gamma * (x_r - centroid)
                std::array<Scalar, Dim> x_e;
                for (size_t d = 0; d < Dim; ++d) {
                    x_e[d] = centroid[d] + gamma * (r.x[d] - centroid[d]);
                }
                Point e = eval(x_e);
                s[Dim] = (e.f < r.f) ? e : r;
            } else if (r.f < s[Dim - 1].f) {
                // Accept reflection
                s[Dim] = r;
            } else {
                // Contraction
                bool outside = (r.f < s[Dim].f);
                std::array<Scalar, Dim> x_c;
                if (outside) {
                    for (size_t d = 0; d < Dim; ++d) {
                        x_c[d] = centroid[d] + rho * (r.x[d] - centroid[d]);
                    }
                } else {
                    for (size_t d = 0; d < Dim; ++d) {
                        x_c[d] = centroid[d] + rho * (s[Dim].x[d] - centroid[d]);
                    }
                }
                Point c = eval(x_c);

                Scalar ref_f = outside ? r.f : s[Dim].f;
                if (c.f < ref_f) {
                    s[Dim] = c;
                } else {
                    // Shrink toward s[0]
                    for (size_t i = 1; i <= Dim; ++i) {
                        std::array<Scalar, Dim> x_s;
                        for (size_t d = 0; d < Dim; ++d) {
                            x_s[d] = s[0].x[d] + sigma * (s[i].x[d] - s[0].x[d]);
                        }
                        s[i] = eval(x_s);
                    }
                }
            }
        }

        std::sort(s.begin(), s.end(), [](const Point& p1, const Point& p2) {
            return p1.f < p2.f;
        });
        return s[0];
    }
};

} // namespace hyphy::opt

