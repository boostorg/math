//  Copyright 2026 Nicholas Thompson and Matt Borland
//  Distributed under the Boost Software License, Version 1.0.
//  https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MATH_QUADRATURE_VAN_DEN_BOS_UNIT_SQUARE_HPP
#define BOOST_MATH_QUADRATURE_VAN_DEN_BOS_UNIT_SQUARE_HPP

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
#error "The header <boost/math/quadrature/van_den_bos_unit_square.hpp> requires C++17 or later."
#endif

#include <boost/math/quadrature/detail/van_den_bos_unit_square_tables.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef BOOST_MATH_INSTRUMENT
#include <iostream>
#endif

namespace boost { namespace math { namespace quadrature {

namespace detail {

struct van_den_bos_unit_square_level_view
{
    std::array<double, 2> const* points;
    double const* weights;
    std::size_t size;
    unsigned degree;
    std::size_t order;
};

inline van_den_bos_unit_square_level_view
van_den_bos_unit_square_level(std::size_t level)
{
    switch (level)
    {
    case 0:
        return {
            van_den_bos_unit_square_level_0_points.data(),
            van_den_bos_unit_square_level_0_weights.data(),
            van_den_bos_unit_square_level_0_points.size(),
            van_den_bos_unit_square_level_0_degree,
            van_den_bos_unit_square_level_0_order};
    case 1:
        return {
            van_den_bos_unit_square_level_1_points.data(),
            van_den_bos_unit_square_level_1_weights.data(),
            van_den_bos_unit_square_level_1_points.size(),
            van_den_bos_unit_square_level_1_degree,
            van_den_bos_unit_square_level_1_order};
    case 2:
        return {
            van_den_bos_unit_square_level_2_points.data(),
            van_den_bos_unit_square_level_2_weights.data(),
            van_den_bos_unit_square_level_2_points.size(),
            van_den_bos_unit_square_level_2_degree,
            van_den_bos_unit_square_level_2_order};
    case 3:
        return {
            van_den_bos_unit_square_level_3_points.data(),
            van_den_bos_unit_square_level_3_weights.data(),
            van_den_bos_unit_square_level_3_points.size(),
            van_den_bos_unit_square_level_3_degree,
            van_den_bos_unit_square_level_3_order};
    case 4:
        return {
            van_den_bos_unit_square_level_4_points.data(),
            van_den_bos_unit_square_level_4_weights.data(),
            van_den_bos_unit_square_level_4_points.size(),
            van_den_bos_unit_square_level_4_degree,
            van_den_bos_unit_square_level_4_order};
    case 5:
        return {
            van_den_bos_unit_square_level_5_points.data(),
            van_den_bos_unit_square_level_5_weights.data(),
            van_den_bos_unit_square_level_5_points.size(),
            van_den_bos_unit_square_level_5_degree,
            van_den_bos_unit_square_level_5_order};
    default:
        throw std::out_of_range(
            "van_den_bos_unit_square: requested level exceeds "
            "the precomputed hierarchy.");
    }
}

static constexpr std::size_t van_den_bos_unit_square_max_level = 5;

} // namespace detail

/*
 * Adaptive positive cubature on the reference square [0,1]^2.
 *
 * The quadrature nodes are a strictly nested immutable hierarchy generated
 * offline in higher precision.  The generator starts from tensor-product
 * Clenshaw-Curtis rules for continuous uniform (Lebesgue) measure and applies
 * protected-node van-den-Bos reduction.  This differs from the equal-weight
 * empirical sample measure used in the paper, while retaining its positive
 * nullspace reduction and nested-node construction.
 *
 * If physical coordinates are required, the caller is responsible for the
 * coordinate map and its Jacobian.
 *
 * Runtime performs no rule construction and mutates no shared state.
 * Independent invocations are therefore thread-safe.
 */
template <class F>
auto van_den_bos_unit_square(
    F const& f,
    double tolerance,
    std::size_t max_level = detail::van_den_bos_unit_square_max_level,
    double* error_estimate = 0,
    double* L1 = 0,
    std::size_t* evaluations = 0)
    -> decltype(f(std::declval<std::array<double, 2> const&>()))
{
    using std::abs;
    typedef std::array<double, 2> point_type;
    typedef decltype(f(std::declval<point_type const&>())) result_type;

    if (!(tolerance > 0))
    {
        throw std::domain_error(
            "van_den_bos_unit_square: tolerance must be positive.");
    }

    max_level = (std::min)(
        max_level, detail::van_den_bos_unit_square_max_level);

    std::vector<result_type> values;

    result_type previous = result_type(0);
    result_type current = result_type(0);

    double current_l1 = 0;
    double current_error =
        (std::numeric_limits<double>::infinity)();

    std::size_t evaluation_count = 0;
    bool have_previous = false;

    for (std::size_t level_index = 0;
         level_index <= max_level;
         ++level_index)
    {
        const detail::van_den_bos_unit_square_level_view level =
            detail::van_den_bos_unit_square_level(level_index);

        // Generator invariant: the previous level is an exact prefix.
        // Therefore only the new suffix requires function evaluations.
        values.reserve(level.size);
        while (values.size() < level.size)
        {
            values.push_back(f(level.points[values.size()]));
            ++evaluation_count;
        }

        current = result_type(0);
        current_l1 = 0;

        for (std::size_t j = 0; j < level.size; ++j)
        {
            current += level.weights[j] * values[j];
            current_l1 += level.weights[j] * abs(values[j]);
        }

#ifdef BOOST_MATH_INSTRUMENT
        std::cout
            << std::setprecision(
                std::numeric_limits<double>::max_digits10)
            << "VDB_LEVEL"
            << ",level=" << level_index
            << ",degree=" << level.degree
            << ",order=" << level.order
            << ",nodes=" << level.size
            << ",evaluations=" << evaluation_count
            << ",value=" << current;
#endif

        if (have_previous)
        {
            const double delta = abs(current - previous);
            const double scale =
                (std::max)(1.0, static_cast<double>(abs(current)));

            // Paper Sec. 2.3: adjacent nested-rule difference.
            // The floor prevents reporting a meaningless zero below binary64
            // roundoff.
            const double roundoff_floor =
                8 * (std::numeric_limits<double>::epsilon)() * scale;
            current_error = (std::max)(delta, roundoff_floor);

#ifdef BOOST_MATH_INSTRUMENT
            std::cout
                << ",delta=" << delta
                << ",error_estimate=" << current_error
                << '\n';
#endif

            if (current_error <= tolerance * scale)
            {
                break;
            }
        }
        else
        {
#ifdef BOOST_MATH_INSTRUMENT
            std::cout << ",delta=nan,error_estimate=nan\n";
#endif
            have_previous = true;
        }

        previous = current;
    }

    if (error_estimate)
    {
        *error_estimate = current_error;
    }
    if (L1)
    {
        *L1 = current_l1;
    }
    if (evaluations)
    {
        *evaluations = evaluation_count;
    }

    return current;
}

}}} // namespaces

#endif
