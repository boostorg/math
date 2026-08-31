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
#include <boost/math/tools/throw_exception.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef BOOST_MATH_INSTRUMENT
#include <iostream>
#include <iomanip>
#endif

namespace boost { namespace math { namespace quadrature {

namespace detail {

struct van_den_bos_unit_square_level_view
{
    std::array<std::uint8_t, 2> const* node_indices;
    double const* weights;
    std::size_t size;

    std::array<double, 2> point(std::size_t index) const
    {
        const std::array<std::uint8_t, 2>& node = node_indices[index];
        return {{
            van_den_bos_unit_square_coordinates[node[0]],
            van_den_bos_unit_square_coordinates[node[1]]}};
    }
};

inline van_den_bos_unit_square_level_view
van_den_bos_unit_square_level(std::size_t level)
{
    if (level >= van_den_bos_unit_square_level_count)
    {
        BOOST_MATH_THROW_EXCEPTION(std::out_of_range(
            "van_den_bos_unit_square: requested level exceeds "
            "the precomputed hierarchy."));
    }

    return {
        van_den_bos_unit_square_node_indices.data(),
        van_den_bos_unit_square_weights[level],
        van_den_bos_unit_square_sizes[level]};
}

} // namespace detail

/*
 * Adaptive positive cubature on the reference square [0,1]^2.
 *
 * The quadrature nodes are a strictly nested immutable hierarchy generated
 * offline in higher precision.  The generator starts from tensor-product
 * endpoint-free Fejer type-II rules for continuous uniform (Lebesgue)
 * measure and applies protected-node van-den-Bos reduction.  This differs
 * from the equal-weight
 * empirical sample measure used in the paper, while retaining its positive
 * nullspace reduction and nested-node construction.
 *
 * For a parallelogram, the caller is responsible for the affine map from
 * [0,1]^2 to physical coordinates and for multiplying by its Jacobian.
 * The rules contain no boundary nodes.
 *
 * Runtime performs no rule construction and mutates no shared state.
 * Independent invocations are therefore thread-safe.
 * Convergence requires two consecutive adjacent-level differences to satisfy
 * the requested tolerance, guarding against an isolated cancellation.
 */
template <class F>
auto van_den_bos_unit_square(
    F const& f,
    double tolerance,
    double* error_estimate = nullptr,
    double* L1 = nullptr,
    std::size_t* evaluations = nullptr)
    -> decltype(f(std::declval<double>(), std::declval<double>()))
{
    using std::abs;
    typedef decltype(
        f(std::declval<double>(), std::declval<double>())) result_type;

    static_assert(!std::is_integral<result_type>::value,
        "The return type cannot be integral, it must be either a real or "
        "complex floating point type.");

    typedef typename std::decay<
        decltype(abs(std::declval<result_type>()))>::type magnitude_type;

    if (!(tolerance > 0))
    {
        BOOST_MATH_THROW_EXCEPTION(std::domain_error(
            "van_den_bos_unit_square: tolerance must be positive."));
    }

    std::vector<result_type> orbit_sums;
    std::vector<magnitude_type> orbit_absolute_sums;

    result_type previous = result_type(0);
    result_type current = result_type(0);

    double current_l1 = 0;
    double current_error =
        (std::numeric_limits<double>::infinity)();

    std::size_t evaluation_count = 0;
    bool have_previous = false;
    bool previous_estimate_satisfied = false;

    for (std::size_t level_index = 0;
         level_index < detail::van_den_bos_unit_square_level_count;
         ++level_index)
    {
        const detail::van_den_bos_unit_square_level_view level =
            detail::van_den_bos_unit_square_level(level_index);

        // Generator invariant: the previous level is an exact prefix.
        // Therefore only the new suffix requires function evaluations.
        orbit_sums.reserve(level.size);
        orbit_absolute_sums.reserve(level.size);
        while (orbit_sums.size() < level.size)
        {
            const std::array<double, 2> point =
                level.point(orbit_sums.size());
            result_type orbit_sum = result_type(0);
            magnitude_type orbit_absolute_sum = magnitude_type(0);
            const auto accumulate = [&](double x, double y)
            {
                const result_type value = f(x, y);
                orbit_sum += value;
                orbit_absolute_sum += abs(value);
                ++evaluation_count;
            };

            const double x = point[0];
            const double y = point[1];
            if (x == 0.5)
            {
                accumulate(x, y);
                if (y != 0.5)
                {
                    accumulate(y, x);
                    accumulate(x, 1 - y);
                    accumulate(1 - y, x);
                }
            }
            else if (x == y)
            {
                accumulate(x, x);
                accumulate(1 - x, x);
                accumulate(x, 1 - x);
                accumulate(1 - x, 1 - x);
            }
            else
            {
                accumulate(x, y);
                accumulate(y, x);
                accumulate(1 - x, y);
                accumulate(y, 1 - x);
                accumulate(x, 1 - y);
                accumulate(1 - y, x);
                accumulate(1 - x, 1 - y);
                accumulate(1 - y, 1 - x);
            }
            orbit_sums.push_back(orbit_sum);
            orbit_absolute_sums.push_back(orbit_absolute_sum);
        }

        current = result_type(0);
        current_l1 = 0;

        for (std::size_t j = 0; j < level.size; ++j)
        {
            current += level.weights[j] * orbit_sums[j];
            current_l1 += level.weights[j] * orbit_absolute_sums[j];
        }

#ifdef BOOST_MATH_INSTRUMENT
        std::cout
            << std::setprecision(
                std::numeric_limits<double>::max_digits10)
            << "VDB_LEVEL"
            << ",level=" << level_index
            << ",nodes=" << level.size
            << ",evaluations=" << evaluation_count
            << ",value=" << current;
#endif

        if (have_previous)
        {
            const double delta = abs(current - previous);
            const double scale =
                (std::max)(1.0, static_cast<double>(abs(current)));

            // Paper Sec. 2.3: adjacent nested-rule difference.  The floor
            // prevents reporting a meaningless zero below the accumulator's
            // precision.
            const double roundoff_floor =
                8 * static_cast<double>(
                    (std::numeric_limits<magnitude_type>::epsilon)()) * scale;
            current_error = (std::max)(delta, roundoff_floor);

#ifdef BOOST_MATH_INSTRUMENT
            std::cout
                << ",delta=" << delta
                << ",error_estimate=" << current_error
                << '\n';
#endif

            const bool estimate_satisfied =
                current_error <= tolerance * scale;

            // One accidentally small adjacent difference is not reliable for
            // non-monotone convergence.  Require confirmation by the next
            // nested rule before stopping.
            if (estimate_satisfied && previous_estimate_satisfied)
            {
                break;
            }
            previous_estimate_satisfied = estimate_satisfied;
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
