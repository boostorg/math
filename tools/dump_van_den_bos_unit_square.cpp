// Copyright 2026 Nicholas Thompson
// Distributed under the Boost Software License, Version 1.0.

/*
 * One-shot generator for nested van den Bos rules on [0,1]^2.
 *
 * Construction precision:
 *   boost::multiprecision::cpp_bin_float_quad
 *
 * Storage/runtime precision:
 *   IEEE-754 binary64, rounded once at emission and written as exact
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
void emit_rule(
    std::ostream& os,
    Rule const& rule,
    std::size_t level,
    unsigned degree,
    std::size_t order)
{
    os << std::hexfloat;
    os
        << "static constexpr unsigned van_den_bos_unit_square_level_" << level
        << "_degree = " << degree << ";\n"
        << "static constexpr std::size_t van_den_bos_unit_square_level_" << level
        << "_order = " << order << ";\n\n";

    os
        << "static constexpr std::array<std::array<double, 2>, "
        << rule.size()
        << "> van_den_bos_unit_square_level_" << level << "_points{{\n";

    for (std::size_t j = 0; j < rule.size(); ++j)
    {
        const double x = static_cast<double>(rule.point(j)[0]);
        const double y = static_cast<double>(rule.point(j)[1]);

        os << "    {{" << x << ", " << y << "}}";
        if (j + 1 != rule.size())
        {
            os << ",";
        }
        os << "\n";
    }
    os << "}};\n\n";

    os
        << "static constexpr std::array<double, "
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

} // namespace

int main()
{
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

    std::cout <<
R"(/*
 * Copyright 2026 Nicholas Thompson and Matt Borland
 * Distributed under the Boost Software License, Version 1.0.
 *
 * Generated by tools/dump_van_den_bos_unit_square.cpp.
 *
 * The Clenshaw-Curtis candidates, moment matrices, support reductions, nodes,
 * and weights were all computed in cpp_bin_float_quad.  Nodes and weights
 * were rounded once, at emission, to IEEE-754 binary64.
 *
 * Every hexadecimal floating-point literal exactly represents its intended
 * IEEE-754 binary64 value.  C++17 or later is therefore required.
 */
#ifndef BOOST_MATH_QUADRATURE_DETAIL_VAN_DEN_BOS_UNIT_SQUARE_TABLES_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_VAN_DEN_BOS_UNIT_SQUARE_TABLES_HPP

#include <array>
#include <cstddef>

namespace boost { namespace math { namespace quadrature { namespace detail {

)";

    for (std::size_t k = 0; k <= 5; ++k)
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

        // Crucially, validate the binary64 rule we are about to freeze,
        // not merely the higher-precision construction.
        validate_rounded_rule(rule, level->degree, k);

        emit_rule(
            std::cout,
            rule,
            k,
            level->degree,
            level->clenshaw_curtis_order);
    }

    std::cout <<
R"(}}}} // namespaces

#endif
)";

    return 0;
}
