/*
 * Copyright 2026 Nicholas Thompson
 *
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0.
 *
 * Positive cubature reduction in the style of:
 *
 * L. M. M. van den Bos, B. Sanderse, W. A. A. M. Bierbooms,
 * G. J. W. van Bussel,
 * "Generating nested quadrature rules with positive weights based on
 * arbitrary sample sets", SIAM/ASA Journal on Uncertainty Quantification.
 */

#ifndef BOOST_MATH_TOOLS_VAN_DEN_BOS_GENERATOR_HPP
#define BOOST_MATH_TOOLS_VAN_DEN_BOS_GENERATOR_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#ifdef BOOST_MATH_INSTRUMENT
#include <iostream>
#endif
#include <memory>
#include <mutex>
#include <queue>
#include <set>
#include <shared_mutex>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#endif

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>

#include <boost/math/tools/precision.hpp>

namespace boost { namespace math { namespace quadrature {

template <class Point, class Real>
class positive_cubature_rule
{
public:
    using point_type = Point;
    using real_type = Real;

    positive_cubature_rule() = default;

    positive_cubature_rule(std::vector<Point> points, std::vector<Real> weights)
        : points_(std::move(points)), weights_(std::move(weights))
    {
        if (points_.size() != weights_.size())
        {
            throw std::invalid_argument(
                "positive_cubature_rule: points and weights must have the same size.");
        }
    }

    std::size_t size() const noexcept
    {
        return points_.size();
    }

    std::vector<Point> const& points() const noexcept
    {
        return points_;
    }

    std::vector<Real> const& weights() const noexcept
    {
        return weights_;
    }

    Point const& point(std::size_t i) const noexcept
    {
        return points_[i];
    }

    Real weight(std::size_t i) const noexcept
    {
        return weights_[i];
    }

    template <class F>
    auto integrate(F const& f, Real* L1 = nullptr) const
        -> decltype(std::declval<F>()(std::declval<Point const&>()))
    {
        using std::abs;
        using result_type = decltype(f(std::declval<Point const&>()));

        static_assert(!std::is_integral<result_type>::value,
            "The integrand return type must be a real or complex floating-point type.");

        result_type result = result_type(0);
        Real l1 = Real(0);

        for (std::size_t i = 0; i < points_.size(); ++i)
        {
            result_type y = f(points_[i]);
            result += weights_[i] * y;
            l1 += weights_[i] * abs(y);
        }

        if (L1)
        {
            *L1 = l1;
        }
        return result;
    }

private:
    std::vector<Point> points_;
    std::vector<Real> weights_;
};

namespace detail {

template <class Real>
using matrix_type = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

template <class Real>
using vector_type = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

template <class Point, class Real, class Basis>
auto fill_moment_values(
    Point const& point,
    Basis const& basis,
    Real* values,
    std::size_t moment_count,
    int)
    -> decltype(basis.evaluate(point, values, moment_count), void())
{
    basis.evaluate(point, values, moment_count);
}

template <class Point, class Real, class Basis>
void fill_moment_values(
    Point const& point,
    Basis const& basis,
    Real* values,
    std::size_t moment_count,
    long)
{
    for (std::size_t k = 0; k < moment_count; ++k)
    {
        values[k] = basis(point, k);
    }
}

template <class Point, class Real, class Basis>
matrix_type<Real> moment_matrix(
    std::vector<Point> const& points,
    std::vector<std::size_t> const& indices,
    Basis const& basis,
    std::size_t moment_count)
{
    matrix_type<Real> A(
        static_cast<Eigen::Index>(moment_count),
        static_cast<Eigen::Index>(indices.size()));

    for (std::size_t j = 0; j < indices.size(); ++j)
    {
        fill_moment_values<Point, Real, Basis>(
            points[indices[j]],
            basis,
            A.col(static_cast<Eigen::Index>(j)).data(),
            moment_count,
            0);
    }
    return A;
}

template <class Point, class Real, class Basis>
matrix_type<Real> moment_matrix(
    std::vector<Point> const& points,
    Basis const& basis,
    std::size_t moment_count)
{
    std::vector<std::size_t> indices(points.size());
    for (std::size_t i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }
    return moment_matrix<Point, Real, Basis>(points, indices, basis, moment_count);
}

template <class Real>
vector_type<Real> to_eigen(std::vector<Real> const& x)
{
    vector_type<Real> result(static_cast<Eigen::Index>(x.size()));
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        result(static_cast<Eigen::Index>(i)) = x[i];
    }
    return result;
}

template <class Real>
std::vector<Real> from_eigen(vector_type<Real> const& x)
{
    std::vector<Real> result(static_cast<std::size_t>(x.size()));
    for (Eigen::Index i = 0; i < x.size(); ++i)
    {
        result[static_cast<std::size_t>(i)] = x(i);
    }
    return result;
}

template <class Real>
matrix_type<Real> kernel(matrix_type<Real> const& A, Real rank_tolerance)
{
    Eigen::ColPivHouseholderQR<matrix_type<Real>> qr(A);
    qr.setThreshold(rank_tolerance);

    const Eigen::Index rank = qr.rank();
    const Eigen::Index n = A.cols();
    const Eigen::Index nullity = n - rank;

    matrix_type<Real> Z = matrix_type<Real>::Zero(n, nullity);
    if (nullity == 0)
    {
        return Z;
    }

    // A P = Q R, with R = [R11 R12].  In permuted coordinates,
    //
    //     R11 y1 + R12 y2 = 0,
    //
    // so a nullspace basis is
    //
    //     y = [-R11^{-1} R12; I].
    matrix_type<Real> R = qr.matrixR().topLeftCorner(rank, n);
    matrix_type<Real> Y = matrix_type<Real>::Zero(n, nullity);

    if (rank > 0)
    {
        Y.topRows(rank) =
            -R.leftCols(rank)
                 .template triangularView<Eigen::Upper>()
                 .solve(R.rightCols(nullity));
    }

    Y.bottomRows(nullity).setIdentity();
    Z = qr.colsPermutation() * Y;
    return Z;
}

template <class Real>
vector_type<Real> kernel_vector(
    matrix_type<Real> const& A,
    Real rank_tolerance)
{
    Eigen::ColPivHouseholderQR<matrix_type<Real>> qr(A);
    qr.setThreshold(rank_tolerance);

    const Eigen::Index rank = qr.rank();
    const Eigen::Index n = A.cols();
    if (rank == n)
    {
        return vector_type<Real>();
    }

    // Set the first free variable in permuted coordinates to one and solve
    // R11 y1 + r12 = 0.  Reduction needs a single null direction, so avoid
    // materializing the potentially very large complete nullspace.
    vector_type<Real> y = vector_type<Real>::Zero(n);
    y(rank) = Real(1);
    if (rank > 0)
    {
        matrix_type<Real> R = qr.matrixR().topLeftCorner(rank, rank + 1);
        y.head(rank) =
            -R.leftCols(rank)
                 .template triangularView<Eigen::Upper>()
                 .solve(R.col(rank));
    }
    return qr.colsPermutation() * y;
}

template <class Real>
bool solve_square(
    matrix_type<Real> const& A,
    std::vector<Real> const& b,
    std::vector<Real>& x,
    Real rank_tolerance)
{
    if (A.rows() != A.cols() || A.rows() != static_cast<Eigen::Index>(b.size()))
    {
        return false;
    }

    Eigen::ColPivHouseholderQR<matrix_type<Real>> qr(A);
    qr.setThreshold(rank_tolerance);
    if (qr.rank() != A.cols())
    {
        return false;
    }

    vector_type<Real> solution = qr.solve(to_eigen(b));
    x = from_eigen(solution);
    return true;
}

template <class Real>
struct reduction_candidate
{
    bool valid = false;
    Real alpha = Real(0);
    Real protected_margin = Real(0);
    std::size_t boundary_index = 0;
};

template <class Real>
reduction_candidate<Real> choose_boundary(
    std::vector<Real> const& w,
    vector_type<Real> const& c,
    std::vector<unsigned char> const& protected_local,
    Real positivity_tolerance)
{
    const Real inf = (std::numeric_limits<Real>::infinity)();

    Real upper = inf;
    Real lower = -inf;
    std::size_t upper_idx = 0;
    std::size_t lower_idx = 0;
    bool have_upper = false;
    bool have_lower = false;

    for (std::size_t i = 0; i < w.size(); ++i)
    {
        Real ci = c(static_cast<Eigen::Index>(i));
        if (ci > positivity_tolerance)
        {
            Real a = w[i] / ci;
            if (!have_upper || a < upper)
            {
                upper = a;
                upper_idx = i;
                have_upper = true;
            }
        }
        else if (ci < -positivity_tolerance)
        {
            Real a = w[i] / ci;
            if (!have_lower || a > lower)
            {
                lower = a;
                lower_idx = i;
                have_lower = true;
            }
        }
    }

    reduction_candidate<Real> best;

    auto inspect = [&](Real alpha, std::size_t boundary_index)
    {
        if (protected_local[boundary_index])
        {
            return;
        }

        Real margin = inf;
        for (std::size_t i = 0; i < w.size(); ++i)
        {
            Real wi = w[i] - alpha * c(static_cast<Eigen::Index>(i));
            if (wi < -positivity_tolerance)
            {
                return;
            }
            if (protected_local[i])
            {
                margin = (std::min)(margin, wi);
            }
        }

        if (!best.valid || margin > best.protected_margin)
        {
            best.valid = true;
            best.alpha = alpha;
            best.protected_margin = margin;
            best.boundary_index = boundary_index;
        }
    };

    if (have_upper)
    {
        inspect(upper, upper_idx);
    }
    if (have_lower)
    {
        inspect(lower, lower_idx);
    }

    return best;
}

inline std::vector<std::size_t> complement_indices(
    std::size_t n,
    std::vector<std::size_t> const& removed)
{
    std::vector<unsigned char> is_removed(n, 0);
    for (std::size_t i : removed)
    {
        is_removed[i] = 1;
    }

    std::vector<std::size_t> result;
    result.reserve(n - removed.size());
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!is_removed[i])
        {
            result.push_back(i);
        }
    }
    return result;
}

template <class Point, class Real, class Basis>
bool weights_for_support(
    std::vector<Point> const& points,
    std::vector<std::size_t> const& support,
    Basis const& basis,
    std::vector<Real> const& moments,
    std::vector<Real>& weights,
    Real rank_tolerance)
{
    if (support.size() != moments.size())
    {
        return false;
    }

    matrix_type<Real> A = moment_matrix<Point, Real, Basis>(
        points, support, basis, moments.size());
    return solve_square(A, moments, weights, rank_tolerance);
}

template <class Point, class Real, class Basis>
bool initial_basic_removal(
    std::vector<Point> const& points,
    std::vector<Real> weights,
    Basis const& basis,
    std::vector<Real> const& moments,
    std::vector<std::size_t>& removed,
    Real rank_tolerance,
    Real positivity_tolerance)
{
    using std::abs;

    std::vector<std::size_t> active(points.size());
    for (std::size_t i = 0; i < active.size(); ++i)
    {
        active[i] = i;
    }
    removed.clear();

    while (active.size() > moments.size())
    {
        matrix_type<Real> A = moment_matrix<Point, Real, Basis>(
            points, active, basis, moments.size());
        matrix_type<Real> Z = kernel(A, rank_tolerance);
        if (Z.cols() == 0)
        {
            return false;
        }

        vector_type<Real> c = Z.col(0);
        Real scale = c.cwiseAbs().maxCoeff();
        if (scale == Real(0))
        {
            return false;
        }
        c /= scale;

        Real upper = (std::numeric_limits<Real>::infinity)();
        Real lower = -(std::numeric_limits<Real>::infinity)();
        bool have_upper = false;
        bool have_lower = false;

        for (std::size_t j = 0; j < active.size(); ++j)
        {
            Real cj = c(static_cast<Eigen::Index>(j));
            Real wj = weights[active[j]];
            if (cj > positivity_tolerance)
            {
                Real a = wj / cj;
                if (!have_upper || a < upper)
                {
                    upper = a;
                    have_upper = true;
                }
            }
            else if (cj < -positivity_tolerance)
            {
                Real a = wj / cj;
                if (!have_lower || a > lower)
                {
                    lower = a;
                    have_lower = true;
                }
            }
        }

        if (!have_upper && !have_lower)
        {
            return false;
        }

        Real alpha;
        if (!have_lower)
        {
            alpha = upper;
        }
        else if (!have_upper)
        {
            alpha = lower;
        }
        else
        {
            alpha = abs(lower) <= abs(upper) ? lower : upper;
        }

        for (std::size_t j = 0; j < active.size(); ++j)
        {
            weights[active[j]] -= alpha * c(static_cast<Eigen::Index>(j));
        }

        bool deleted = false;
        std::vector<std::size_t> next;
        next.reserve(active.size() - 1);
        for (std::size_t index : active)
        {
            if (!deleted && abs(weights[index]) <= positivity_tolerance)
            {
                weights[index] = Real(0);
                removed.push_back(index);
                deleted = true;
            }
            else
            {
                next.push_back(index);
            }
        }

        if (!deleted)
        {
            return false;
        }
        active.swap(next);
    }

    std::sort(removed.begin(), removed.end());
    return true;
}


template <std::size_t Dimension>
std::vector<std::array<unsigned, Dimension>>
total_degree_multi_indices(unsigned degree)
{
    std::vector<std::array<unsigned, Dimension>> result;
    std::array<unsigned, Dimension> alpha{};

    auto append = [&](auto&& self,
                      std::size_t coordinate,
                      unsigned remaining) -> void
    {
        if (coordinate + 1 == Dimension)
        {
            alpha[coordinate] = remaining;
            result.push_back(alpha);
            return;
        }

        for (unsigned k = 0; k <= remaining; ++k)
        {
            alpha[coordinate] = k;
            self(self, coordinate + 1, remaining - k);
        }
    };

    for (unsigned total = 0; total <= degree; ++total)
    {
        append(append, 0, total);
    }
    return result;
}

template <class Real>
Real shifted_legendre(unsigned n, Real x)
{
    const Real t = Real(2) * x - Real(1);

    if (n == 0)
    {
        return Real(1);
    }
    if (n == 1)
    {
        return t;
    }

    Real pnm2 = Real(1);
    Real pnm1 = t;
    for (unsigned k = 2; k <= n; ++k)
    {
        const Real pn =
            (Real(2 * k - 1) * t * pnm1 - Real(k - 1) * pnm2) / Real(k);
        pnm2 = pnm1;
        pnm1 = pn;
    }
    return pnm1;
}

/*
 * Reduce only the unprotected suffix to a basic positive solution while the
 * protected prefix and its current candidate weights remain fixed.  This is
 * the straightforward nested extension described at the start of paper
 * Section 3.2, before the minimal M-removal construction.
 *
 * A basis inverse is reused across simplex-style pivots.  Rebuilding it
 * periodically avoids the cost of a fresh rectangular QR decomposition for
 * every removed node while keeping the offline multiprecision construction
 * numerically controlled.
 */
template <class Point, class Real, class Basis>
positive_cubature_rule<Point, Real> reduce_protected_basic(
    std::vector<Point> points,
    std::vector<Real> weights,
    Basis const& basis,
    std::vector<Real> const& moments,
    std::size_t protected_count,
    Real rank_tolerance,
    Real positivity_tolerance)
{
    using std::abs;

    const std::size_t moment_count = moments.size();
    const std::size_t suffix_count = points.size() - protected_count;
    if (suffix_count < moment_count)
    {
        throw std::runtime_error(
            "reduce_protected_basic: too few removable candidate nodes.");
    }

    std::vector<std::size_t> suffix_indices(suffix_count);
    for (std::size_t j = 0; j < suffix_count; ++j)
    {
        suffix_indices[j] = protected_count + j;
    }

    matrix_type<Real> A = moment_matrix<Point, Real, Basis>(
        points, suffix_indices, basis, moment_count);

    std::vector<Real> residual = moments;
    std::vector<Real> protected_moments(moment_count);
    for (std::size_t j = 0; j < protected_count; ++j)
    {
        fill_moment_values<Point, Real, Basis>(
            points[j],
            basis,
            protected_moments.data(),
            moment_count,
            0);
        for (std::size_t k = 0; k < moment_count; ++k)
        {
            residual[k] -= weights[j] * protected_moments[k];
        }
    }

    Eigen::ColPivHouseholderQR<matrix_type<Real>> qr(A);
    qr.setThreshold(rank_tolerance);
    if (qr.rank() != static_cast<Eigen::Index>(moment_count))
    {
        throw std::runtime_error(
            "reduce_protected_basic: candidate moment matrix is rank deficient.");
    }

    std::vector<std::size_t> basic(moment_count);
    std::vector<unsigned char> is_basic(suffix_count, 0);
    for (std::size_t j = 0; j < moment_count; ++j)
    {
        const std::size_t index = static_cast<std::size_t>(
            qr.colsPermutation().indices()(static_cast<Eigen::Index>(j)));
        basic[j] = index;
        is_basic[index] = 1;
    }

    matrix_type<Real> B(
        static_cast<Eigen::Index>(moment_count),
        static_cast<Eigen::Index>(moment_count));
    matrix_type<Real> inverse;
    const Real solve_tolerance = Real(1000) * rank_tolerance;

    auto inverse_residual = [&]
    {
        return (B * inverse - matrix_type<Real>::Identity(
                    static_cast<Eigen::Index>(moment_count),
                    static_cast<Eigen::Index>(moment_count)))
            .cwiseAbs()
            .maxCoeff();
    };

    auto rebuild_inverse = [&]
    {
        for (std::size_t j = 0; j < moment_count; ++j)
        {
            B.col(static_cast<Eigen::Index>(j)) =
                A.col(static_cast<Eigen::Index>(basic[j]));
        }

        Eigen::PartialPivLU<matrix_type<Real>> basis_lu(B);
        inverse = basis_lu.inverse();

        if (inverse_residual() > solve_tolerance)
        {
            Eigen::ColPivHouseholderQR<matrix_type<Real>> basis_qr(B);
            basis_qr.setThreshold(rank_tolerance);
            if (basis_qr.rank() != static_cast<Eigen::Index>(moment_count))
            {
                throw std::runtime_error(
                    "reduce_protected_basic: pivot basis is rank deficient.");
            }
            inverse = basis_qr.solve(matrix_type<Real>::Identity(
                static_cast<Eigen::Index>(moment_count),
                static_cast<Eigen::Index>(moment_count)));
            if (inverse_residual() > solve_tolerance)
            {
                throw std::runtime_error(
                    "reduce_protected_basic: inverse residual is too large.");
            }
        }
    };
    rebuild_inverse();

    std::vector<Real> basic_weights(moment_count);
    std::vector<Real> nonbasic_weights(suffix_count);
    for (std::size_t j = 0; j < suffix_count; ++j)
    {
        nonbasic_weights[j] = weights[protected_count + j];
    }
    for (std::size_t j = 0; j < moment_count; ++j)
    {
        basic_weights[j] = weights[protected_count + basic[j]];
        nonbasic_weights[basic[j]] = Real(0);
    }

    for (std::size_t candidate = 0; candidate < suffix_count; ++candidate)
    {
        if (is_basic[candidate])
        {
            continue;
        }

        if (abs(nonbasic_weights[candidate]) <= positivity_tolerance)
        {
            nonbasic_weights[candidate] = Real(0);
            continue;
        }

        const auto candidate_column =
            A.col(static_cast<Eigen::Index>(candidate));
        vector_type<Real> z = inverse * candidate_column;
        if ((B * z - candidate_column).cwiseAbs().maxCoeff()
            > solve_tolerance)
        {
            rebuild_inverse();
            z = inverse * candidate_column;
            if ((B * z - candidate_column).cwiseAbs().maxCoeff()
                > solve_tolerance)
            {
                throw std::runtime_error(
                    "reduce_protected_basic: candidate solve residual is too large.");
            }
        }
        const Real candidate_weight = nonbasic_weights[candidate];

        struct boundary
        {
            bool valid = false;
            Real alpha = Real(0);
            std::size_t killed = 0;
            std::size_t zero_count = 0;
        };

        auto inspect = [&](Real alpha, std::size_t killed)
        {
            boundary result;
            result.valid = true;
            result.alpha = alpha;
            result.killed = killed;

            const Real updated_candidate = candidate_weight - alpha;
            if (abs(updated_candidate) <= positivity_tolerance)
            {
                ++result.zero_count;
            }
            else if (updated_candidate < -positivity_tolerance)
            {
                result.valid = false;
            }

            for (std::size_t j = 0; j < moment_count; ++j)
            {
                const Real updated =
                    basic_weights[j] +
                    alpha * z(static_cast<Eigen::Index>(j));
                if (abs(updated) <= positivity_tolerance)
                {
                    ++result.zero_count;
                }
                else if (updated < -positivity_tolerance)
                {
                    result.valid = false;
                }
            }
            return result;
        };

        Real upper = candidate_weight;
        std::size_t upper_killed = moment_count;
        Real lower = -(std::numeric_limits<Real>::infinity)();
        std::size_t lower_killed = moment_count;
        bool have_lower = false;

        for (std::size_t j = 0; j < moment_count; ++j)
        {
            const Real zj = z(static_cast<Eigen::Index>(j));
            if (zj < -positivity_tolerance)
            {
                const Real ratio = basic_weights[j] / (-zj);
                if (ratio < upper)
                {
                    upper = ratio;
                    upper_killed = j;
                }
            }
            else if (zj > positivity_tolerance)
            {
                const Real ratio = -basic_weights[j] / zj;
                if (!have_lower || ratio > lower)
                {
                    lower = ratio;
                    lower_killed = j;
                    have_lower = true;
                }
            }
        }

        boundary chosen = inspect(upper, upper_killed);
        if (have_lower)
        {
            boundary alternative = inspect(lower, lower_killed);
            if (!chosen.valid ||
                alternative.zero_count < chosen.zero_count ||
                (alternative.zero_count == chosen.zero_count &&
                 abs(alternative.alpha) < abs(chosen.alpha)))
            {
                chosen = alternative;
            }
        }

        if (!chosen.valid || chosen.zero_count != 1)
        {
            throw std::runtime_error(
                "reduce_protected_basic: degenerate positivity boundary.");
        }

        std::vector<Real> updated_weights(moment_count);
        for (std::size_t j = 0; j < moment_count; ++j)
        {
            updated_weights[j] =
                basic_weights[j] +
                chosen.alpha * z(static_cast<Eigen::Index>(j));
            if (abs(updated_weights[j]) <= positivity_tolerance)
            {
                updated_weights[j] = Real(0);
            }
        }
        Real updated_candidate = candidate_weight - chosen.alpha;
        if (abs(updated_candidate) <= positivity_tolerance)
        {
            updated_candidate = Real(0);
        }

        if (chosen.killed == moment_count)
        {
            // The candidate left the active support; the basis is unchanged.
            basic_weights = std::move(updated_weights);
            nonbasic_weights[candidate] = Real(0);
            continue;
        }

        const std::size_t replaced = chosen.killed;
        if (!(updated_candidate > Real(0)) ||
            z(static_cast<Eigen::Index>(replaced)) == Real(0))
        {
            throw std::runtime_error(
                "reduce_protected_basic: invalid basis pivot.");
        }

        const std::size_t removed = basic[replaced];
        is_basic[removed] = 0;
        nonbasic_weights[removed] = Real(0);
        is_basic[candidate] = 1;
        nonbasic_weights[candidate] = Real(0);
        basic[replaced] = candidate;
        B.col(static_cast<Eigen::Index>(replaced)) = candidate_column;
        basic_weights = std::move(updated_weights);
        basic_weights[replaced] = updated_candidate;

        vector_type<Real> update = z;
        update(static_cast<Eigen::Index>(replaced)) -= Real(1);
        inverse -=
            (update * inverse.row(static_cast<Eigen::Index>(replaced))) /
            z(static_cast<Eigen::Index>(replaced));

    }

    std::sort(basic.begin(), basic.end());
    std::vector<std::size_t> support(moment_count);
    for (std::size_t j = 0; j < moment_count; ++j)
    {
        support[j] = protected_count + basic[j];
    }

    std::vector<Real> final_weights;
    if (!weights_for_support(
            points,
            support,
            basis,
            residual,
            final_weights,
            rank_tolerance))
    {
        throw std::runtime_error(
            "reduce_protected_basic: final support is singular.");
    }

    std::vector<Point> reduced_points;
    std::vector<Real> reduced_weights;
    reduced_points.reserve(protected_count + moment_count);
    reduced_weights.reserve(protected_count + moment_count);

    for (std::size_t j = 0; j < protected_count; ++j)
    {
        reduced_points.push_back(points[j]);
        reduced_weights.push_back(weights[j]);
    }
    for (std::size_t j = 0; j < moment_count; ++j)
    {
        if (!(final_weights[j] > Real(0)))
        {
            throw std::runtime_error(
                "reduce_protected_basic: final weight is not positive.");
        }
        reduced_points.push_back(points[support[j]]);
        reduced_weights.push_back(final_weights[j]);
    }

    return positive_cubature_rule<Point, Real>(
        std::move(reduced_points), std::move(reduced_weights));
}

} // namespace detail

template <class Point, class Real, class Basis>
positive_cubature_rule<Point, Real> van_den_bos_reduce(
    std::vector<Point> points,
    std::vector<Real> weights,
    Basis const& basis,
    std::vector<Real> const& moments,
    std::size_t target_size = 0,
    std::size_t protected_count = 0,
    Real rank_tolerance = Real(100) * boost::math::tools::epsilon<Real>(),
    Real positivity_tolerance = Real(1000) * boost::math::tools::epsilon<Real>())
{
    using std::abs;

    if (points.size() != weights.size())
    {
        throw std::invalid_argument(
            "van_den_bos_reduce: points and weights must have the same size.");
    }
    if (protected_count > points.size())
    {
        throw std::invalid_argument(
            "van_den_bos_reduce: protected_count exceeds number of points.");
    }
    if (moments.empty())
    {
        throw std::invalid_argument(
            "van_den_bos_reduce: moments must not be empty.");
    }

    for (Real w : weights)
    {
        if (w < -positivity_tolerance)
        {
            throw std::domain_error(
                "van_den_bos_reduce: input rule must have nonnegative weights.");
        }
    }

    if (target_size == 0)
    {
        target_size = moments.size();
    }
    target_size = (std::max)(target_size, protected_count);

    std::vector<unsigned char> protected_flag(points.size(), 0);
    for (std::size_t i = 0; i < protected_count; ++i)
    {
        protected_flag[i] = 1;
    }

    for (std::size_t i = points.size(); i-- > protected_count;)
    {
        if (abs(weights[i]) <= positivity_tolerance)
        {
            points.erase(points.begin() + static_cast<std::ptrdiff_t>(i));
            weights.erase(weights.begin() + static_cast<std::ptrdiff_t>(i));
            protected_flag.erase(
                protected_flag.begin() + static_cast<std::ptrdiff_t>(i));
        }
    }

    if (protected_count > 0 &&
        target_size == protected_count + moments.size())
    {
        return detail::reduce_protected_basic(
            std::move(points),
            std::move(weights),
            basis,
            moments,
            protected_count,
            rank_tolerance,
            positivity_tolerance);
    }

    while (points.size() > target_size)
    {
        // Compute directions supported only on removable nodes.  Consequently
        // the protected prefix retains its coordinates and positive weights
        // while each boundary step deletes one unprotected node and preserves
        // every requested moment.
        std::vector<std::size_t> removable(points.size() - protected_count);
        for (std::size_t i = 0; i < removable.size(); ++i)
        {
            removable[i] = protected_count + i;
        }

        detail::matrix_type<Real> A =
            detail::moment_matrix<Point, Real, Basis>(
                points, removable, basis, moments.size());
        detail::vector_type<Real> local_c =
            detail::kernel_vector(A, rank_tolerance);
        if (local_c.size() == 0)
        {
            break;
        }

        bool reduced = false;
        Real scale = local_c.cwiseAbs().maxCoeff();
        if (scale != Real(0))
        {
            local_c /= scale;

            detail::vector_type<Real> c =
                detail::vector_type<Real>::Zero(
                    static_cast<Eigen::Index>(points.size()));
            for (std::size_t i = 0; i < removable.size(); ++i)
            {
                c(static_cast<Eigen::Index>(removable[i])) =
                    local_c(static_cast<Eigen::Index>(i));
            }

            auto candidate = detail::choose_boundary(
                weights, c, protected_flag, positivity_tolerance);
            if (candidate.valid)
            {
                for (std::size_t i = 0; i < weights.size(); ++i)
                {
                    weights[i] -=
                        candidate.alpha * c(static_cast<Eigen::Index>(i));
                    if (abs(weights[i]) <= positivity_tolerance)
                    {
                        weights[i] = Real(0);
                    }
                }

                if (protected_count == 0)
                {
                    // An unprotected basic reduction may hit a higher
                    // dimensional face.  Remove every boundary node, as in
                    // Algorithm 1; the exact-size postcondition below rejects
                    // an unsuitable reduction path.
                    for (std::size_t i = points.size(); i-- > 0;)
                    {
                        if (abs(weights[i]) <= positivity_tolerance)
                        {
                            points.erase(
                                points.begin() +
                                static_cast<std::ptrdiff_t>(i));
                            weights.erase(
                                weights.begin() +
                                static_cast<std::ptrdiff_t>(i));
                            protected_flag.erase(
                                protected_flag.begin() +
                                static_cast<std::ptrdiff_t>(i));
                        }
                    }
                }
                else
                {
                    const std::size_t removed = candidate.boundary_index;
                    if (protected_flag[removed])
                    {
                        throw std::runtime_error(
                            "van_den_bos_reduce: attempted to remove a protected node.");
                    }
                    points.erase(
                        points.begin() + static_cast<std::ptrdiff_t>(removed));
                    weights.erase(
                        weights.begin() + static_cast<std::ptrdiff_t>(removed));
                    protected_flag.erase(
                        protected_flag.begin() +
                        static_cast<std::ptrdiff_t>(removed));
                }
                reduced = true;
            }
        }

        if (!reduced)
        {
            break;
        }
    }

    if (points.size() != target_size)
    {
        throw std::runtime_error(
            "van_den_bos_reduce: could not reach the requested nested support size.");
    }

    return positive_cubature_rule<Point, Real>(
        std::move(points), std::move(weights));
}

/*
 * Enumerate the M-removals of van den Bos et al.  If n nodes and r moment
 * constraints are supplied, M=n-r.  Each returned sorted set contains M
 * nodes whose simultaneous removal leaves a positive r-node cubature rule.
 */
template <class Point, class Real, class Basis>
std::vector<std::vector<std::size_t>> van_den_bos_m_removals(
    std::vector<Point> const& points,
    std::vector<Real> const& weights,
    Basis const& basis,
    std::vector<Real> const& moments,
    std::size_t max_removals = 0,
    Real rank_tolerance = Real(100) * boost::math::tools::epsilon<Real>(),
    Real positivity_tolerance = Real(1000) * boost::math::tools::epsilon<Real>())
{
    using std::abs;

    if (points.size() != weights.size() || moments.empty() ||
        points.size() < moments.size())
    {
        throw std::invalid_argument("van_den_bos_m_removals: invalid sizes.");
    }

    const std::size_t n = points.size();
    const std::size_t r = moments.size();
    const std::size_t M = n - r;

    if (M == 0)
    {
        return {{}};
    }

    std::vector<std::size_t> initial;
    if (!detail::initial_basic_removal(
            points, weights, basis, moments, initial,
            rank_tolerance, positivity_tolerance))
    {
        return {};
    }

    std::queue<std::vector<std::size_t>> queue;
    std::set<std::vector<std::size_t>> seen;
    std::vector<std::vector<std::size_t>> result;

    queue.push(initial);
    seen.insert(initial);

    while (!queue.empty())
    {
        std::vector<std::size_t> q = queue.front();
        queue.pop();
        result.push_back(q);

        if (max_removals && result.size() >= max_removals)
        {
            break;
        }

        for (std::size_t ri = 0; ri < M; ++ri)
        {
            std::vector<std::size_t> fixed_removed;
            fixed_removed.reserve(M - 1);
            for (std::size_t j = 0; j < M; ++j)
            {
                if (j != ri)
                {
                    fixed_removed.push_back(q[j]);
                }
            }

            std::vector<std::size_t> extended_support =
                detail::complement_indices(n, fixed_removed);
            detail::matrix_type<Real> A =
                detail::moment_matrix<Point, Real, Basis>(
                    points, extended_support, basis, r);
            detail::matrix_type<Real> Z = detail::kernel(A, rank_tolerance);
            if (Z.cols() != 1)
            {
                continue;
            }
            detail::vector_type<Real> c = Z.col(0);

            std::vector<std::size_t> support =
                detail::complement_indices(n, q);
            std::vector<Real> basic_weights;
            if (!detail::weights_for_support(
                    points, support, basis, moments,
                    basic_weights, rank_tolerance))
            {
                continue;
            }

            std::vector<Real> extended_weights(extended_support.size(), Real(0));
            std::size_t b = 0;
            std::size_t restored_local = extended_support.size();
            for (std::size_t j = 0; j < extended_support.size(); ++j)
            {
                if (extended_support[j] == q[ri])
                {
                    restored_local = j;
                }
                else
                {
                    extended_weights[j] = basic_weights[b++];
                }
            }
            if (restored_local == extended_support.size())
            {
                continue;
            }

            Real upper = (std::numeric_limits<Real>::infinity)();
            Real lower = -(std::numeric_limits<Real>::infinity)();
            bool have_upper = false;
            bool have_lower = false;
            std::size_t upper_index = 0;
            std::size_t lower_index = 0;

            for (std::size_t j = 0; j < extended_support.size(); ++j)
            {
                Real cj = c(static_cast<Eigen::Index>(j));
                if (cj > positivity_tolerance)
                {
                    Real a = extended_weights[j] / cj;
                    if (!have_upper || a < upper)
                    {
                        upper = a;
                        upper_index = j;
                        have_upper = true;
                    }
                }
                else if (cj < -positivity_tolerance)
                {
                    Real a = extended_weights[j] / cj;
                    if (!have_lower || a > lower)
                    {
                        lower = a;
                        lower_index = j;
                        have_lower = true;
                    }
                }
            }

            std::size_t killed = extended_support.size();
            if (have_upper && abs(upper) > positivity_tolerance)
            {
                killed = upper_index;
            }
            if (have_lower && abs(lower) > positivity_tolerance)
            {
                if (killed == extended_support.size() || abs(lower) > abs(upper))
                {
                    killed = lower_index;
                }
            }
            if (killed == extended_support.size())
            {
                continue;
            }

            std::size_t replacement = extended_support[killed];
            if (replacement == q[ri])
            {
                continue;
            }

            std::vector<std::size_t> qhat = q;
            qhat[ri] = replacement;
            std::sort(qhat.begin(), qhat.end());
            if (std::adjacent_find(qhat.begin(), qhat.end()) != qhat.end())
            {
                continue;
            }

            std::vector<std::size_t> candidate_support =
                detail::complement_indices(n, qhat);
            std::vector<Real> candidate_weights;
            if (!detail::weights_for_support(
                    points, candidate_support, basis, moments,
                    candidate_weights, rank_tolerance))
            {
                continue;
            }

            bool positive = true;
            for (Real w : candidate_weights)
            {
                if (w < -positivity_tolerance)
                {
                    positive = false;
                    break;
                }
            }
            if (!positive)
            {
                continue;
            }

            if (seen.insert(qhat).second)
            {
                queue.push(std::move(qhat));
            }
        }
    }

    return result;
}

template <class Point, class Real, class Basis>
positive_cubature_rule<Point, Real> van_den_bos_reduce_multiple(
    std::vector<Point> const& points,
    std::vector<Real> const& weights,
    Basis const& basis,
    std::vector<Real> const& moments,
    std::size_t protected_count = 0,
    bool require_protected = true,
    std::size_t max_removals = 0,
    Real rank_tolerance = Real(100) * boost::math::tools::epsilon<Real>(),
    Real positivity_tolerance = Real(1000) * boost::math::tools::epsilon<Real>())
{
    if (protected_count > points.size())
    {
        throw std::invalid_argument(
            "van_den_bos_reduce_multiple: invalid protected_count.");
    }

    auto removals = van_den_bos_m_removals(
        points, weights, basis, moments, max_removals,
        rank_tolerance, positivity_tolerance);
    if (removals.empty())
    {
        throw std::runtime_error(
            "van_den_bos_reduce_multiple: no M-removal found.");
    }

    std::vector<std::size_t> const* best = nullptr;
    std::size_t best_hits = (std::numeric_limits<std::size_t>::max)();

    for (auto const& q : removals)
    {
        std::size_t hits = 0;
        for (std::size_t i : q)
        {
            if (i < protected_count)
            {
                ++hits;
            }
        }

        if (require_protected && hits != 0)
        {
            continue;
        }
        if (!best || hits < best_hits || (hits == best_hits && q < *best))
        {
            best = &q;
            best_hits = hits;
        }
    }

    if (!best)
    {
        throw std::runtime_error(
            "van_den_bos_reduce_multiple: no removal preserves protected nodes.");
    }

    std::vector<std::size_t> support =
        detail::complement_indices(points.size(), *best);
    std::vector<Real> final_weights;
    if (!detail::weights_for_support(
            points, support, basis, moments,
            final_weights, rank_tolerance))
    {
        throw std::runtime_error(
            "van_den_bos_reduce_multiple: selected support is singular.");
    }

    std::vector<Point> final_points;
    final_points.reserve(support.size());
    for (std::size_t i : support)
    {
        final_points.push_back(points[i]);
    }

    return positive_cubature_rule<Point, Real>(
        std::move(final_points), std::move(final_weights));
}

template <class Point, class Real, class Basis>
Real van_den_bos_moment_residual(
    positive_cubature_rule<Point, Real> const& rule,
    Basis const& basis,
    std::vector<Real> const& moments)
{
    using std::abs;

    Real result = 0;
    for (std::size_t k = 0; k < moments.size(); ++k)
    {
        Real q = 0;
        for (std::size_t j = 0; j < rule.size(); ++j)
        {
            q += rule.weight(j) * basis(rule.point(j), k);
        }
        result = (std::max)(result, abs(q - moments[k]));
    }
    return result;
}


/*
 * Total-degree tensor-product shifted Legendre basis on [0,1]^Dimension.
 *
 * van den Bos et al. use products of Legendre polynomials rather than
 * monomials to improve conditioning of the Vandermonde matrix.  For
 * Lebesgue measure on the unit cube, every nonconstant basis function has
 * zero moment, so the exact moment vector is (1, 0, ..., 0).
 */
template <class Real, std::size_t Dimension>
class van_den_bos_legendre_basis
{
public:
    using point_type = std::array<Real, Dimension>;
    using multi_index_type = std::array<unsigned, Dimension>;

    explicit van_den_bos_legendre_basis(unsigned degree)
        : degree_(degree),
          indices_(detail::total_degree_multi_indices<Dimension>(degree))
    {}

    std::size_t size() const noexcept
    {
        return indices_.size();
    }

    Real operator()(point_type const& x, std::size_t k) const
    {
        Real value = Real(1);
        auto const& alpha = indices_[k];
        for (std::size_t j = 0; j < Dimension; ++j)
        {
            value *= detail::shifted_legendre(alpha[j], x[j]);
        }
        return value;
    }

    void evaluate(
        point_type const& x,
        Real* values,
        std::size_t count) const
    {
        std::array<std::vector<Real>, Dimension> polynomials;
        for (std::size_t j = 0; j < Dimension; ++j)
        {
            polynomials[j].resize(degree_ + 1);
            polynomials[j][0] = Real(1);
            if (degree_ == 0)
            {
                continue;
            }

            const Real t = Real(2) * x[j] - Real(1);
            polynomials[j][1] = t;
            for (unsigned k = 2; k <= degree_; ++k)
            {
                polynomials[j][k] =
                    (Real(2 * k - 1) * t * polynomials[j][k - 1]
                     - Real(k - 1) * polynomials[j][k - 2]) / Real(k);
            }
        }

        for (std::size_t k = 0; k < count; ++k)
        {
            Real value = Real(1);
            for (std::size_t j = 0; j < Dimension; ++j)
            {
                value *= polynomials[j][indices_[k][j]];
            }
            values[k] = value;
        }
    }

    std::vector<Real> moments() const
    {
        std::vector<Real> result(size(), Real(0));
        if (!result.empty())
        {
            result[0] = Real(1);
        }
        return result;
    }

    std::vector<multi_index_type> const& indices() const noexcept
    {
        return indices_;
    }

private:
    unsigned degree_;
    std::vector<multi_index_type> indices_;
};

// Legendre products invariant under the full symmetry group of the square.
// Orbit weights are per-node weights, so each basis value includes the size
// of the node's D4 orbit.
template <class Real>
class van_den_bos_d4_legendre_basis
{
public:
    using point_type = std::array<Real, 2>;
    using index_type = std::array<unsigned, 2>;

    explicit van_den_bos_d4_legendre_basis(unsigned degree)
    {
        for (unsigned a = 0; a <= degree; a += 2)
        {
            for (unsigned b = a; b <= degree - a; b += 2)
            {
                indices_.push_back({{a, b}});
            }
        }
    }

    std::size_t size() const noexcept
    {
        return indices_.size();
    }

    Real operator()(point_type const& point, std::size_t k) const
    {
        if (k == 0)
        {
            return static_cast<Real>(orbit_size(point));
        }

        const auto& index = indices_[k];
        const auto transformed = transforms(point);
        Real sum = Real(0);
        for (auto const& p : transformed)
        {
            sum += detail::shifted_legendre(index[0], p[0]) *
                detail::shifted_legendre(index[1], p[1]);
        }
        return sum * static_cast<Real>(orbit_size(point)) /
            static_cast<Real>(8);
    }

    std::vector<Real> moments() const
    {
        std::vector<Real> result(size(), Real(0));
        if (!result.empty())
        {
            result[0] = Real(1);
        }
        return result;
    }

private:
    static std::vector<point_type> transforms(point_type const& p)
    {
        return {
            {{p[0], p[1]}}, {{p[1], p[0]}},
            {{Real(1) - p[0], p[1]}},
            {{p[1], Real(1) - p[0]}},
            {{p[0], Real(1) - p[1]}},
            {{Real(1) - p[1], p[0]}},
            {{Real(1) - p[0], Real(1) - p[1]}},
            {{Real(1) - p[1], Real(1) - p[0]}}};
    }

    static std::size_t orbit_size(point_type const& p)
    {
        const auto transformed = transforms(p);
        std::size_t size = 0;
        for (std::size_t i = 0; i < transformed.size(); ++i)
        {
            bool duplicate = false;
            for (std::size_t j = 0; j < i; ++j)
            {
                bool equal = true;
                for (std::size_t d = 0; d < 2; ++d)
                {
                    Real scale = (std::max)(Real(1), (std::max)(
                        abs(transformed[i][d]), abs(transformed[j][d])));
                    if (abs(transformed[i][d] - transformed[j][d]) >
                        Real(64) * std::numeric_limits<Real>::epsilon() * scale)
                        equal = false;
                }
                if (equal)
                {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate)
            {
                ++size;
            }
        }
        return size;
    }

    std::vector<index_type> indices_;
};


namespace detail {

template <class Real>
std::pair<std::vector<Real>, std::vector<Real>>
clenshaw_curtis_rule(std::size_t n)
{
    using std::acos;
    using std::cos;

    if (n == 0)
    {
        return {{Real(0.5)}, {Real(1)}};
    }

    const Real pi = acos(Real(-1));
    std::vector<Real> x(n + 1);
    std::vector<Real> w(n + 1);

    for (std::size_t j = 0; j <= n; ++j)
    {
        const Real theta = pi * Real(j) / Real(n);
        // Map the usual [-1,1] Clenshaw-Curtis nodes to [0,1].
        x[j] = (Real(1) - cos(theta)) / Real(2);
    }

    if (n == 1)
    {
        w[0] = Real(0.5);
        w[1] = Real(0.5);
        return {x, w};
    }

    if ((n & 1U) == 0U)
    {
        const Real endpoint = Real(1) / (Real(2) * Real(n * n - 1));
        w[0] = endpoint;
        w[n] = endpoint;

        for (std::size_t j = 1; j < n; ++j)
        {
            const Real theta = pi * Real(j) / Real(n);
            Real v = Real(1);

            for (std::size_t k = 1; k < n / 2; ++k)
            {
                v -= Real(2) * cos(Real(2 * k) * theta)
                    / Real(4 * k * k - 1);
            }

            v -= cos(Real(n) * theta) / Real(n * n - 1);
            w[j] = v / Real(n);
        }
    }
    else
    {
        const Real endpoint = Real(1) / (Real(2) * Real(n * n));
        w[0] = endpoint;
        w[n] = endpoint;

        for (std::size_t j = 1; j < n; ++j)
        {
            const Real theta = pi * Real(j) / Real(n);
            Real v = Real(1);

            for (std::size_t k = 1; k <= (n - 1) / 2; ++k)
            {
                v -= Real(2) * cos(Real(2 * k) * theta)
                    / Real(4 * k * k - 1);
            }

            w[j] = v / Real(n);
        }
    }

    return {x, w};
}

template <class Real, std::size_t Dimension>
std::pair<
    std::vector<std::array<Real, Dimension>>,
    std::vector<Real>>
clenshaw_curtis_tensor_rule(std::size_t n)
{
    using point_type = std::array<Real, Dimension>;

    auto one_dimensional = clenshaw_curtis_rule<Real>(n);
    auto const& x = one_dimensional.first;
    auto const& w = one_dimensional.second;

    std::size_t count = 1;
    for (std::size_t d = 0; d < Dimension; ++d)
    {
        count *= n + 1;
    }

    std::vector<point_type> points;
    std::vector<Real> weights;
    points.reserve(count);
    weights.reserve(count);

    std::array<std::size_t, Dimension> index{};
    bool done = false;

    while (!done)
    {
        point_type point{};
        Real weight = Real(1);

        for (std::size_t d = 0; d < Dimension; ++d)
        {
            point[d] = x[index[d]];
            weight *= w[index[d]];
        }

        points.push_back(point);
        weights.push_back(weight);

        for (std::size_t d = Dimension; d-- > 0;)
        {
            ++index[d];
            if (index[d] <= n)
            {
                break;
            }

            index[d] = 0;
            if (d == 0)
            {
                done = true;
            }
        }
    }

    return {std::move(points), std::move(weights)};
}

// Nested Fejer type-II nodes.  The parameter is the denominator in
// theta_j = j*pi/order, so powers of two form a nested sequence while the
// endpoints are excluded.
template <class Real>
std::pair<std::vector<Real>, std::vector<Real>>
fejer_second_rule(std::size_t order)
{
    using std::acos;
    using std::cos;
    using std::sin;

    if (order < 2)
        throw std::invalid_argument("fejer_second_rule: order must be at least 2");

    const Real pi = acos(Real(-1));
    std::vector<Real> x(order - 1);
    std::vector<Real> w(order - 1);
    for (std::size_t j = 1; j < order; ++j)
    {
        const Real theta = pi * Real(j) / Real(order);
        x[j - 1] = (Real(1) - cos(theta)) / Real(2);
        Real sum = Real(0);
        for (std::size_t k = 1; k <= order / 2; ++k)
        {
            const std::size_t odd = 2 * k - 1;
            sum += sin(Real(odd) * theta) / Real(odd);
        }
        // The factor is two on [-1,1] and one after mapping to [0,1].
        w[j - 1] = Real(2) * sin(theta) * sum / Real(order);
    }
    return {std::move(x), std::move(w)};
}

template <class Real, std::size_t Dimension>
std::pair<std::vector<std::array<Real, Dimension>>, std::vector<Real>>
fejer_second_tensor_rule(std::size_t order)
{
    using point_type = std::array<Real, Dimension>;
    auto one_dimensional = fejer_second_rule<Real>(order);
    auto const& x = one_dimensional.first;
    auto const& w = one_dimensional.second;
    std::size_t count = 1;
    for (std::size_t d = 0; d < Dimension; ++d) count *= x.size();
    std::vector<point_type> points;
    std::vector<Real> weights;
    points.reserve(count);
    weights.reserve(count);
    std::array<std::size_t, Dimension> index{};
    bool done = false;
    while (!done)
    {
        point_type point{};
        Real weight = Real(1);
        for (std::size_t d = 0; d < Dimension; ++d)
        {
            point[d] = x[index[d]];
            weight *= w[index[d]];
        }
        points.push_back(point);
        weights.push_back(weight);
        for (std::size_t d = Dimension; d-- > 0;)
        {
            if (++index[d] < x.size()) break;
            index[d] = 0;
            if (d == 0) done = true;
        }
    }
    return {std::move(points), std::move(weights)};
}

template <class Real, std::size_t Dimension>
bool same_point(
    std::array<Real, Dimension> const& x,
    std::array<Real, Dimension> const& y)
{
    using std::abs;

    const Real tolerance =
        Real(64) * (std::numeric_limits<Real>::epsilon)();

    for (std::size_t d = 0; d < Dimension; ++d)
    {
        if (abs(x[d] - y[d]) >
            tolerance * (std::max)(Real(1), (std::max)(abs(x[d]), abs(y[d]))))
        {
            return false;
        }
    }
    return true;
}

inline std::size_t next_power_of_two(std::size_t n)
{
    std::size_t result = 1;
    while (result < n)
    {
        if (result > (std::numeric_limits<std::size_t>::max)() / 2)
        {
            throw std::overflow_error(
                "van_den_bos_cache: candidate-grid order overflow.");
        }
        result *= 2;
    }
    return result;
}

template <class Real, std::size_t Dimension>
void place_protected_prefix(
    positive_cubature_rule<std::array<Real, Dimension>, Real> const& old_rule,
    std::vector<std::array<Real, Dimension>>& points,
    std::vector<Real>& weights)
{
    using point_type = std::array<Real, Dimension>;

    std::vector<point_type> reordered_points;
    std::vector<Real> reordered_weights;
    reordered_points.reserve(points.size());
    reordered_weights.reserve(weights.size());

    std::vector<unsigned char> used(points.size(), 0);

    for (point_type const& old_point : old_rule.points())
    {
        std::size_t match = points.size();

        for (std::size_t j = 0; j < points.size(); ++j)
        {
            if (!used[j] && same_point(old_point, points[j]))
            {
                match = j;
                break;
            }
        }

        if (match == points.size())
        {
            throw std::runtime_error(
                "van_den_bos_cache: refined candidate grid does not contain an old node.");
        }

        used[match] = 1;
        // Preserve the old point bit-for-bit so cached evaluations can be
        // carried forward without coordinate comparisons.
        reordered_points.push_back(old_point);
        reordered_weights.push_back(weights[match]);
    }

    for (std::size_t j = 0; j < points.size(); ++j)
    {
        if (!used[j])
        {
            reordered_points.push_back(points[j]);
            reordered_weights.push_back(weights[j]);
        }
    }

    points = std::move(reordered_points);
    weights = std::move(reordered_weights);
}

} // namespace detail

/*
 * Thread-safe shared hierarchy of immutable van-den-Bos cubature rules.
 *
 * Individual integrations remain serial.  Many integrations may share one
 * cache concurrently.  Reading an existing level takes a shared lock; only
 * construction of a new level takes the exclusive lock.  The second check
 * under the exclusive lock prevents duplicate refinement when several threads
 * request the same new level simultaneously.
 *
 * Refinement uses dyadically nested Clenshaw-Curtis tensor candidates on
 * [0,1]^Dimension and the shifted-Legendre moment basis recommended by
 * van den Bos et al.  Polynomial-space dimension is increased by at most
 * max_moment_growth at each level.  Each reduced rule contains the preceding
 * rule as an exact prefix, followed by at most one new node per moment.
 */
template <class Real, std::size_t Dimension>
class van_den_bos_cache
{
public:
    using real_type = Real;
    using point_type = std::array<Real, Dimension>;
    using rule_type = positive_cubature_rule<point_type, Real>;

    struct level_type
    {
        std::shared_ptr<rule_type const> rule;
        unsigned degree = 0;
        std::size_t clenshaw_curtis_order = 0;
    };

    explicit van_den_bos_cache(
        unsigned initial_degree = 4,
        std::size_t max_moment_growth = 2,
        Real rank_tolerance =
            Real(100) * boost::math::tools::epsilon<Real>(),
        Real positivity_tolerance =
            Real(1000) * boost::math::tools::epsilon<Real>())
        : initial_degree_(initial_degree),
          max_moment_growth_(max_moment_growth),
          rank_tolerance_(rank_tolerance),
          positivity_tolerance_(positivity_tolerance)
    {
        if (initial_degree_ == 0)
        {
            throw std::invalid_argument(
                "van_den_bos_cache: initial_degree must be positive.");
        }
        if (max_moment_growth_ < 2)
        {
            throw std::invalid_argument(
                "van_den_bos_cache: max_moment_growth must be at least 2.");
        }
        levels_.push_back(make_level(nullptr, initial_degree_));
    }

    van_den_bos_cache(
        rule_type initial_rule,
        unsigned initial_degree,
        std::size_t max_moment_growth = 2,
        Real rank_tolerance =
            Real(100) * boost::math::tools::epsilon<Real>(),
        Real positivity_tolerance =
            Real(1000) * boost::math::tools::epsilon<Real>())
        : initial_degree_(initial_degree),
          max_moment_growth_(max_moment_growth),
          rank_tolerance_(rank_tolerance),
          positivity_tolerance_(positivity_tolerance)
    {
        if (initial_rule.size() == 0 || initial_degree_ == 0)
        {
            throw std::invalid_argument(
                "van_den_bos_cache: initial rule and degree must be nonzero.");
        }
        if (max_moment_growth_ < 2)
        {
            throw std::invalid_argument(
                "van_den_bos_cache: max_moment_growth must be at least 2.");
        }

        auto initial = std::make_shared<level_type>();
        initial->rule = std::make_shared<rule_type const>(
            std::move(initial_rule));
        initial->degree = initial_degree_;
        initial->clenshaw_curtis_order =
            detail::next_power_of_two(
                static_cast<std::size_t>(initial_degree_));
        levels_.push_back(initial);
    }

    std::shared_ptr<level_type const> level(std::size_t index) const
    {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        if (index >= levels_.size())
        {
            return {};
        }
        return levels_[index];
    }

    std::shared_ptr<level_type const> ensure_level(std::size_t index)
    {
        {
            std::shared_lock<std::shared_mutex> lock(mutex_);
            if (index < levels_.size())
            {
                return levels_[index];
            }
        }

        std::unique_lock<std::shared_mutex> lock(mutex_);
        while (index >= levels_.size())
        {
            auto const& previous = levels_.back();
            const unsigned next_degree = choose_next_degree(*previous);
            levels_.push_back(make_level(previous.get(), next_degree));
        }
        return levels_[index];
    }

    std::size_t size() const
    {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        return levels_.size();
    }

private:
    unsigned choose_next_degree(level_type const& previous) const
    {
        van_den_bos_legendre_basis<Real, Dimension> previous_basis(
            previous.degree);
        const std::size_t budget =
            max_moment_growth_ * previous_basis.size();

        unsigned degree = previous.degree + 1;
        while (true)
        {
            van_den_bos_legendre_basis<Real, Dimension> candidate(degree + 1);
            if (candidate.size() > budget)
            {
                break;
            }
            ++degree;
        }
        return degree;
    }

    std::shared_ptr<level_type const> make_level(
        level_type const* previous,
        unsigned degree) const
    {
        van_den_bos_legendre_basis<Real, Dimension> basis(degree);
        const std::size_t order = detail::next_power_of_two(
            static_cast<std::size_t>(degree));
        auto moments = basis.moments();

        auto candidate =
            detail::clenshaw_curtis_tensor_rule<Real, Dimension>(order);
        auto points = std::move(candidate.first);
        auto weights = std::move(candidate.second);

        const std::size_t protected_count =
            previous ? previous->rule->size() : 0;
        if (previous)
        {
            detail::place_protected_prefix(
                *previous->rule, points, weights);
        }

        const std::size_t target_size = protected_count + basis.size();
        rule_type reduced = van_den_bos_reduce(
            std::move(points),
            std::move(weights),
            basis,
            moments,
            target_size,
            protected_count,
            rank_tolerance_,
            positivity_tolerance_);

        if (previous)
        {
            for (std::size_t i = 0; i < protected_count; ++i)
            {
                if (reduced.point(i) != previous->rule->point(i))
                {
                    throw std::runtime_error(
                        "van_den_bos_cache: nested prefix invariant failed.");
                }
            }
        }

        auto result = std::make_shared<level_type>();
        result->rule =
            std::make_shared<rule_type const>(std::move(reduced));
        result->degree = degree;
        result->clenshaw_curtis_order = order;
        return result;
    }

    unsigned initial_degree_;
    std::size_t max_moment_growth_;
    Real rank_tolerance_;
    Real positivity_tolerance_;

    mutable std::shared_mutex mutex_;
    std::vector<std::shared_ptr<level_type const>> levels_;
};

/*
 * Adaptive integration over the unit cube [0,1]^Dimension.
 *
 * One call is deliberately serial.  The shared cache only synchronizes rule
 * refinement; function values are private to this invocation.  If another
 * thread has already published the next level, this call immediately reuses
 * that immutable rule.
 *
 * Function values are cached by point inside this invocation.  A later rule
 * may have a different support, but any point already sampled is reused.
 * This avoids a global integrand-value cache, which would be incorrect when
 * many stiffness-matrix entries use different integrands.
 */
template <class F, class Real, std::size_t Dimension>
auto van_den_bos(
    F const& f,
    van_den_bos_cache<Real, Dimension>& cache,
    Real tolerance,
    std::size_t max_refinements = 10,
    Real* error_estimate = nullptr,
    Real* L1 = nullptr,
    std::size_t* evaluations = nullptr)
    -> decltype(f(std::declval<std::array<Real, Dimension> const&>()))
{
    using std::abs;
    using result_type =
        decltype(f(std::declval<std::array<Real, Dimension> const&>()));

    if (!(tolerance > Real(0)))
    {
        throw std::domain_error(
            "van_den_bos: tolerance must be positive.");
    }

    std::vector<std::array<Real, Dimension>> sampled_points;
    std::vector<result_type> sampled_values;
    result_type previous = result_type(0);
    result_type current = result_type(0);
    Real current_l1 = Real(0);
    Real current_error = (std::numeric_limits<Real>::infinity)();
    Real previous_delta = (std::numeric_limits<Real>::infinity)();
    std::size_t evaluation_count = 0;
    std::shared_ptr<typename van_den_bos_cache<Real, Dimension>::rule_type const>
        previous_rule;

    auto same_rule = [](auto const& lhs, auto const& rhs)
    {
        using std::abs;

        if (lhs.size() != rhs.size())
        {
            return false;
        }

        const Real tolerance =
            Real(256) * (std::numeric_limits<Real>::epsilon)();

        for (std::size_t j = 0; j < lhs.size(); ++j)
        {
            if (!detail::same_point(lhs.point(j), rhs.point(j)))
            {
                return false;
            }

            const Real scale =
                (std::max)(Real(1),
                    (std::max)(abs(lhs.weight(j)), abs(rhs.weight(j))));
            if (abs(lhs.weight(j) - rhs.weight(j)) > tolerance * scale)
            {
                return false;
            }
        }
        return true;
    };

    for (std::size_t level_index = 0;
         level_index <= max_refinements;
         ++level_index)
    {
        auto level = cache.ensure_level(level_index);
        auto const& rule = *level->rule;

        current = result_type(0);
        current_l1 = Real(0);

        for (std::size_t j = 0; j < rule.size(); ++j)
        {
            auto const& point = rule.point(j);

            std::size_t cached = sampled_points.size();
            for (std::size_t k = 0; k < sampled_points.size(); ++k)
            {
                if (detail::same_point(point, sampled_points[k]))
                {
                    cached = k;
                    break;
                }
            }

            if (cached == sampled_points.size())
            {
                sampled_points.push_back(point);
                sampled_values.push_back(f(point));
                ++evaluation_count;
            }

            auto const& value = sampled_values[cached];
            current += rule.weight(j) * value;
            current_l1 += rule.weight(j) * abs(value);
        }

#ifdef BOOST_MATH_INSTRUMENT
        std::cout
            << "VDB_LEVEL"
            << ",level=" << level_index
            << ",degree=" << level->degree
            << ",order=" << level->clenshaw_curtis_order
            << ",nodes=" << rule.size()
            << ",evaluations=" << evaluation_count
            << ",value=" << current;
#endif

        if (previous_rule && !same_rule(*previous_rule, rule))
        {
            const Real delta = abs(current - previous);
            const Real scale =
                (std::max)(Real(1), static_cast<Real>(abs(current)));

            // With three genuinely different rules, estimate the remaining
            // geometric tail from the observed contraction:
            //
            //   r_k = delta_k / delta_{k-1}
            //   error_k ~= r_k * delta_k / (1 - r_k)
            //
            // If contraction is absent, fall back to the adjacent-rule
            // difference.
            current_error = delta;
            if (std::isfinite(previous_delta) &&
                previous_delta > Real(0))
            {
                const Real ratio = delta / previous_delta;
                if (ratio >= Real(0) && ratio < Real(1))
                {
                    current_error =
                        ratio * delta / (Real(1) - ratio);
                }
            }

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

            previous_delta = delta;
            previous = current;
            previous_rule = level->rule;
        }
        else if (!previous_rule)
        {
#ifdef BOOST_MATH_INSTRUMENT
            std::cout
                << ",delta=nan"
                << ",error_estimate=nan"
                << '\n';
#endif
            previous = current;
            previous_rule = level->rule;
        }
        else
        {
#ifdef BOOST_MATH_INSTRUMENT
            std::cout
                << ",delta=0"
                << ",error_estimate=" << current_error
                << '\n';
#endif
        }
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
