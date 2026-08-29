// Copyright 2026 Nicholas Thompson
// Distributed under the Boost Software License, Version 1.0.

/*
 * One-shot generator for nested van den Bos rules on [0,1]^2.
 *
 * Construction precision:
 *   boost::multiprecision::cpp_bin_float_quad
 *
 * Storage/runtime precision:
 *   double, rounded once at emission and written as exact
 *   hexadecimal floating-point literals.
 *
 * Build from the Boost tree:
 *
 * c++ -O2 -DNDEBUG -std=c++17 -Iinclude -I/path/to/eigen \
 *   tools/dump_van_den_bos_unit_square.cpp -o dump_vdb_unit_square
 *
 * ./dump_vdb_unit_square \
 *   > include/boost/math/quadrature/detail/van_den_bos_unit_square_tables.hpp
 */

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/eigen.hpp>

#include "van_den_bos_generator.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using quad = boost::multiprecision::cpp_bin_float_quad;
using point_quad = std::array<quad, 2>;
using point_double = std::array<double, 2>;

template <class Rule>
void validate_quad_rule(Rule const& rule, unsigned degree, std::size_t level)
{
    using boost::math::quadrature::van_den_bos_legendre_basis;
    using boost::math::quadrature::van_den_bos_moment_residual;
    using std::abs;

    van_den_bos_legendre_basis<quad, 2> basis(degree);
    auto moments = basis.moments();

    quad sum_weights = 0;
    quad min_weight = (std::numeric_limits<quad>::infinity)();
    for (std::size_t j = 0; j < rule.size(); ++j)
    {
        const quad w = rule.weight(j);
        if (!(w > 0))
        {
            throw std::runtime_error(
                "nonpositive quad weight at level " +
                std::to_string(level));
        }
        sum_weights += w;
        min_weight = (std::min)(min_weight, w);

        for (std::size_t d = 0; d < 2; ++d)
        {
            const quad x = rule.point(j)[d];
            if (x < 0 || x > 1)
            {
                throw std::runtime_error(
                    "quad node outside the unit square at level " +
                    std::to_string(level));
            }
        }
    }

    const quad sum_error = abs(sum_weights - 1);
    const quad residual =
        van_den_bos_moment_residual(rule, basis, moments);
    const quad tolerance =
        quad(1000000000) * std::numeric_limits<quad>::epsilon();

    std::cerr
        << "level " << level
        << ": quad sumw error=" << std::scientific
        << std::setprecision(std::numeric_limits<quad>::max_digits10)
        << sum_error
        << ", quad minw=" << min_weight
        << ", quad moment residual=" << residual
        << '\n';

    if (sum_error > tolerance || residual > tolerance)
    {
        throw std::runtime_error(
            "quad rule moment residual is too large at level " +
            std::to_string(level));
    }
}

template <class Rule>
double rounded_moment_residual(Rule const& rule, unsigned degree)
{
    using boost::math::quadrature::van_den_bos_legendre_basis;

    van_den_bos_legendre_basis<double, 2> basis(degree);
    auto moments = basis.moments();

    double max_residual = 0;

    for (std::size_t k = 0; k < basis.size(); ++k)
    {
        double moment = 0;
        for (std::size_t j = 0; j < rule.size(); ++j)
        {
            point_double p{
                static_cast<double>(rule.point(j)[0]),
                static_cast<double>(rule.point(j)[1])
            };
            const double w = static_cast<double>(rule.weight(j));
            moment += w * basis(p, k);
        }

        max_residual =
            (std::max)(max_residual, std::abs(moment - moments[k]));
    }

    return max_residual;
}

template <class Rule>
void validate_rounded_rule(Rule const& rule, unsigned degree, std::size_t level)
{
    double sum_weights = 0;
    double min_weight = (std::numeric_limits<double>::infinity)();

    for (std::size_t j = 0; j < rule.size(); ++j)
    {
        const double w = static_cast<double>(rule.weight(j));

        if (!std::isfinite(w))
        {
            throw std::runtime_error(
                "non-finite rounded weight at level " + std::to_string(level));
        }

        sum_weights += w;
        min_weight = (std::min)(min_weight, w);

        for (std::size_t d = 0; d < 2; ++d)
        {
            const double x = static_cast<double>(rule.point(j)[d]);
            if (!std::isfinite(x) || x < 0.0 || x > 1.0)
            {
                throw std::runtime_error(
                    "invalid rounded node at level " + std::to_string(level));
            }
        }

        for (std::size_t i = 0; i < j; ++i)
        {
            if (static_cast<double>(rule.point(i)[0]) ==
                    static_cast<double>(rule.point(j)[0]) &&
                static_cast<double>(rule.point(i)[1]) ==
                    static_cast<double>(rule.point(j)[1]))
            {
                throw std::runtime_error(
                    "duplicate rounded node at level " +
                    std::to_string(level));
            }
        }
    }

    const double residual = rounded_moment_residual(rule, degree);

    std::cerr
        << "level " << level
        << ": rounded sumw error=" << std::setprecision(17)
        << std::abs(sum_weights - 1.0)
        << ", rounded minw=" << min_weight
        << ", rounded moment residual=" << residual
        << '\n';

    if (min_weight < 0.0)
    {
        throw std::runtime_error(
            "negative rounded weight at level " + std::to_string(level));
    }

    if (std::abs(sum_weights - 1.0) > 1e-13)
    {
        throw std::runtime_error(
            "rounded weights do not sum to one at level "
            + std::to_string(level));
    }

    if (residual > 5e-13)
    {
        throw std::runtime_error(
            "rounded rule moment residual is too large at level "
            + std::to_string(level));
    }
}

template <class Rule>
void emit_level(
    std::ostream& os,
    Rule const& rule,
    std::size_t level)
{
    os << std::hexfloat;
    os
        << "inline constexpr std::array<double, "
        << rule.size()
        << "> van_den_bos_unit_square_level_" << level << "_weights{{\n";

    for (std::size_t j = 0; j < rule.size(); ++j)
    {
        const double w = static_cast<double>(rule.weight(j));

        os << "    " << w;
        if (j + 1 != rule.size())
        {
            os << ",";
        }
        os << "\n";
    }
    os << "}};\n\n";
}

template <class Rule>
void emit_points(std::ostream& os, Rule const& rule)
{
    os << std::hexfloat;
    std::vector<double> coordinates;
    std::vector<std::array<std::size_t, 2>> indices;
    indices.reserve(rule.size());

    for (std::size_t j = 0; j < rule.size(); ++j)
    {
        std::array<std::size_t, 2> point_indices{};
        for (std::size_t d = 0; d < 2; ++d)
        {
            const double coordinate =
                static_cast<double>(rule.point(j)[d]);
            auto const found =
                std::find(coordinates.begin(), coordinates.end(), coordinate);
            if (found == coordinates.end())
            {
                point_indices[d] = coordinates.size();
                coordinates.push_back(coordinate);
            }
            else
            {
                point_indices[d] = static_cast<std::size_t>(
                    found - coordinates.begin());
            }
        }
        indices.push_back(point_indices);
    }

    if (coordinates.size() >
        static_cast<std::size_t>((std::numeric_limits<std::uint8_t>::max)()) + 1)
    {
        throw std::runtime_error(
            "too many distinct coordinates for uint8_t node indices");
    }

    os
        << "inline constexpr std::array<double, " << coordinates.size()
        << "> van_den_bos_unit_square_coordinates{{\n";

    for (std::size_t j = 0; j < coordinates.size(); ++j)
    {
        os << "    " << coordinates[j]
           << (j + 1 == coordinates.size() ? "\n" : ",\n");
    }
    os << "}};\n\n";

    os
        << "inline constexpr std::array<std::array<std::uint8_t, 2>, "
        << rule.size()
        << "> van_den_bos_unit_square_node_indices{{\n";

    for (std::size_t j = 0; j < rule.size(); ++j)
    {
        os << "    {{" << indices[j][0] << ", " << indices[j][1] << "}}";
        if (j + 1 != rule.size())
        {
            os << ",";
        }
        os << "\n";
    }
    os << "}};\n\n";
}

void emit_hierarchy(std::ostream& os, std::size_t max_level)
{
    const std::size_t level_count = max_level + 1;
    os
        << "inline constexpr std::size_t "
        << "van_den_bos_unit_square_level_count = " << level_count
        << ";\n\n";

    os << "inline constexpr std::array<double const*, " << level_count
       << "> van_den_bos_unit_square_weights{{\n";
    for (std::size_t k = 0; k < level_count; ++k)
    {
        os << "    van_den_bos_unit_square_level_" << k << "_weights.data()"
           << (k + 1 == level_count ? "\n" : ",\n");
    }
    os << "}};\n\n";

    os << "inline constexpr std::array<std::size_t, " << level_count
       << "> van_den_bos_unit_square_sizes{{\n";
    for (std::size_t k = 0; k < level_count; ++k)
    {
        os << "    van_den_bos_unit_square_level_" << k << "_weights.size()"
           << (k + 1 == level_count ? "\n" : ",\n");
    }
    os << "}};\n\n";
}

void emit_header(std::ostream& os)
{
    os <<
R"(/*
 * Copyright 2026 Nicholas Thompson and Matt Borland
 * Distributed under the Boost Software License, Version 1.0.
 *
 * Generated by tools/dump_van_den_bos_unit_square.cpp.
 *
 * The Fejer type-II candidates, moment matrices, support reductions, nodes,
 * and weights were all computed in cpp_bin_float_quad.  Nodes and weights
 * were rounded once, at emission, to double.
 *
 * Every hexadecimal floating-point literal exactly represents its intended
 * double value.  C++17 or later is therefore required.
 */
#ifndef BOOST_MATH_QUADRATURE_DETAIL_VAN_DEN_BOS_UNIT_SQUARE_TABLES_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_VAN_DEN_BOS_UNIT_SQUARE_TABLES_HPP

#include <array>
#include <cstddef>
#include <cstdint>

namespace boost { namespace math { namespace quadrature { namespace detail {

)";
}

void emit_footer(std::ostream& os)
{
    os <<
R"(}}}} // namespaces

#endif
)";
}

} // namespace

// Experimental symmetry-constrained search; this does not alter emitted tables.
int run_symmetric(std::size_t max_level)
{
    using namespace boost::math::quadrature;
    using rule_type = positive_cubature_rule<point_quad, quad>;
    std::vector<rule_type> levels;
    std::vector<rule_type> full_levels;
    emit_header(std::cout);
    unsigned degree = 4;
    for (std::size_t k = 0; k <= max_level; ++k)
    {
        van_den_bos_d4_legendre_basis<quad> basis(degree);
        const std::size_t protected_count = k ? levels.back().size() : 0;
        const std::size_t target_size = protected_count + basis.size();
        std::size_t order = detail::next_power_of_two(degree + 1);
        std::vector<point_quad> points;
        std::vector<quad> weights;
        // Keep only fundamental-triangle representatives.  If this grid is
        // too small for the protected prefix plus the new moment space,
        // advance to the next nested Fejer grid.
        while (true)
        {
            auto candidate = detail::fejer_second_tensor_rule<quad, 2>(order);
            points.clear();
            weights.clear();
            for (std::size_t i = 0; i < candidate.first.size(); ++i)
            {
                auto p = candidate.first[i];
                if (p[1] > p[0] || p[0] > quad("0.5")) continue;
                points.push_back(p);
                weights.push_back(candidate.second[i]);
            }
            if (points.size() >= target_size) break;
            order *= 2;
        }
        if (k) detail::place_protected_prefix(levels.back(), points, weights);
        auto reduced = van_den_bos_reduce(std::move(points), std::move(weights),
            basis, basis.moments(), target_size, protected_count,
            quad(100) * boost::math::tools::epsilon<quad>(),
            quad(1000) * boost::math::tools::epsilon<quad>());
        std::vector<point_quad> full_points;
        std::vector<quad> full_weights;
        for (std::size_t i = 0; i < reduced.size(); ++i)
        {
            auto p = reduced.point(i);
            std::array<point_quad, 8> orbit{{
                {{p[0],p[1]}}, {{p[1],p[0]}}, {{quad(1)-p[0],p[1]}},
                {{p[1],quad(1)-p[0]}}, {{p[0],quad(1)-p[1]}},
                {{quad(1)-p[1],p[0]}}, {{quad(1)-p[0],quad(1)-p[1]}},
                {{quad(1)-p[1],quad(1)-p[0]}}}};
            for (auto const& q : orbit)
            {
                bool duplicate = false;
                for (auto const& old : full_points)
                    if (detail::same_point(q, old)) { duplicate = true; break; }
                if (!duplicate) { full_points.push_back(q); full_weights.push_back(reduced.weight(i)); }
            }
        }
        rule_type full(std::move(full_points), std::move(full_weights));
        for (std::size_t i = 0; i < full.size(); ++i)
            for (std::size_t d = 0; d < 2; ++d)
            {
                const quad x = full.point(i)[d];
                const double rounded_x = static_cast<double>(x);
                if (!(x > 0 && x < 1) ||
                    !(rounded_x > 0 && rounded_x < 1))
                    throw std::runtime_error(
                        "endpoint node in symmetric Fejer construction");
            }
        validate_quad_rule(full, degree, k);
        validate_rounded_rule(full, degree, k);
        if (k)
        {
            if (full.size() <= full_levels.back().size())
                throw std::runtime_error("D4 hierarchy did not grow");
            for (std::size_t i = 0; i < full_levels.back().size(); ++i)
                if (!detail::same_point(full.point(i), full_levels.back().point(i)))
                    throw std::runtime_error("D4 full-node prefix failed");
        }
        levels.push_back(std::move(reduced));
        full_levels.push_back(std::move(full));
        emit_level(std::cout, levels.back(), k);
        std::cerr << "symmetric level " << k << ": degree " << degree
                  << ", representatives " << levels.back().size()
                  << ", full nodes " << full_levels.back().size() << "\n";
        std::size_t budget = 2 * basis.size();
        ++degree;
        while (van_den_bos_d4_legendre_basis<quad>(degree + 1).size() <= budget) ++degree;
    }
    emit_points(std::cout, levels.back());
    emit_hierarchy(std::cout, max_level);
    emit_footer(std::cout);
    return 0;
}

int main(int argc, char** argv)
{
    if (argc == 1)
        return run_symmetric(6);
    if (std::string(argv[1]) == "--symmetric")
    {
        std::size_t max_level = argc >= 3 ? static_cast<std::size_t>(std::stoul(argv[2])) : 6;
        return run_symmetric(max_level);
    }
    if (std::string(argv[1]) != "--legacy")
        return run_symmetric(static_cast<std::size_t>(std::stoul(argv[1])));
    --argc;
    ++argv;
    using boost::math::quadrature::positive_cubature_rule;
    using boost::math::quadrature::van_den_bos_cache;

    // The full construction is performed in cpp_bin_float_quad.  A full
    // tensor Clenshaw-Curtis rule is the reproducible, strictly positive seed.
    // It is already exact for the initial total-degree-four space, so no
    // numerically fragile minimal-support seed is needed.
    auto seed =
        boost::math::quadrature::detail::
            clenshaw_curtis_tensor_rule<quad, 2>(4);

    van_den_bos_cache<quad, 2> cache(
        positive_cubature_rule<point_quad, quad>(
            std::move(seed.first), std::move(seed.second)),
        4);

    emit_header(std::cout);

    std::size_t max_level = 6;
    if (argc == 2)
    {
        max_level = static_cast<std::size_t>(std::stoul(argv[1]));
        if (max_level > 6)
        {
            throw std::invalid_argument("maximum generated level is 6");
        }
    }
    else if (argc != 1)
    {
        throw std::invalid_argument("usage: dump_vdb_unit_square [max-level]");
    }

    for (std::size_t k = 0; k <= max_level; ++k)
    {
        auto level = cache.ensure_level(k);
        auto const& rule = *level->rule;

        if (k > 0)
        {
            auto previous = cache.ensure_level(k - 1);
            if (rule.size() <= previous->rule->size())
            {
                throw std::runtime_error(
                    "nested hierarchy did not grow at level " +
                    std::to_string(k));
            }

            for (std::size_t j = 0; j < previous->rule->size(); ++j)
            {
                for (std::size_t d = 0; d < 2; ++d)
                {
                    if (rule.point(j)[d] != previous->rule->point(j)[d])
                    {
                        throw std::runtime_error(
                            "quad prefix nesting failed at level " +
                            std::to_string(k));
                    }
                    if (static_cast<double>(rule.point(j)[d]) !=
                        static_cast<double>(previous->rule->point(j)[d]))
                    {
                        throw std::runtime_error(
                            "rounded prefix nesting failed at level " +
                            std::to_string(k));
                    }
                }
            }
        }

        std::cerr
            << "constructed level " << k
            << ": degree=" << level->degree
            << ", order=" << level->clenshaw_curtis_order
            << ", nodes=" << rule.size()
            << '\n';

        validate_quad_rule(rule, level->degree, k);

        // Crucially, validate the double rule we are about to freeze,
        // not merely the higher-precision construction.
        validate_rounded_rule(rule, level->degree, k);

        emit_level(std::cout, rule, k);
    }

    emit_points(std::cout, *cache.ensure_level(max_level)->rule);
    emit_hierarchy(std::cout, max_level);

    emit_footer(std::cout);

    return 0;
}
