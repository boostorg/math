// Copyright 2014 Marco Guazzone (marco.guazzone@gmail.com).
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <algorithm>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/hyperexponential.hpp>
#include <boost/math/tools/precision.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cstddef>
#include <iostream>
#include <vector>

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
typedef boost::mpl::list<float, double, long double, boost::math::concepts::real_concept> test_types;
#else
typedef boost::mpl::list<float, double> test_types;
#endif

template <typename RealT>
RealT make_tolerance()
{
    // We must use a low precision since it seems the involved computations are very challenging from the numerical point of view.
    // Indeed, both Octave 3.6.4, MATLAB 2012a and Mathematica 10 provides different results.
    // E.g.:
    //  x = [0 1 2 3 4]
    //  p = [0.2 0.3 0.5]
    //  r = [0.5 1.0 1.5]
    //  PDF(x)
    //    - MATLAB:      1.033333333333333,  0.335636985323608,  0.135792553231720,   0.061039382459897,   0.028790027125382
    //    - Octave:      1.0333333333333332, 0.3356369853236084, 0.1357925532317197,  0.0610393824598966,  0.0287900271253818
    //    - Mathematica: 1.15,               0.3383645184340184, 0.11472883036402601, 0.04558088392888389, 0.02088728412278129
    //
    //  (Tested under Fedora Linux 20 x86_64 running on Intel(R) Core(TM) i7-3540M)
    //

/*
    RealT tol = std::max(boost::math::tools::epsilon<RealT>(),
                         static_cast<RealT>(boost::math::tools::epsilon<double>()*5)*150);

    // At float precision we need to up the tolerance, since 
    // the input values are rounded off to inexact quantities
    // the results get thrown off by a noticeable amount.

    if (boost::math::tools::digits<RealT>() < 50)
    {
        tol *= 50;
    }
    if (boost::is_floating_point<RealT>::value != 1)
    {
        tol *= 20; // real_concept special functions are less accurate
    }
*/
   // Current test data is limited to double precision:
   const RealT tol = (std::max)(static_cast<RealT>(boost::math::tools::epsilon<double>()), boost::math::tools::epsilon<RealT>()) * 100 * 100;

    //std::cout << "[" << __func__ << "] Tolerance: " << tol << "%" << std::endl;

    return tol;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(range, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs) / sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    std::pair<RealT,RealT> res;
    res = boost::math::range(dist);

    BOOST_CHECK_CLOSE( res.first, static_cast<RealT>(0), tol );
    if(std::numeric_limits<RealT>::has_infinity)
    {
       BOOST_CHECK_EQUAL(res.second, std::numeric_limits<RealT>::infinity());
    }
    else
    {
       BOOST_CHECK_EQUAL(res.second, boost::math::tools::max_value<RealT>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(support, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs)/sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    std::pair<RealT,RealT> res;
    res = boost::math::support(dist);

    BOOST_CHECK_CLOSE( res.first, boost::math::tools::min_value<RealT>(), tol );
    BOOST_CHECK_CLOSE( res.second, boost::math::tools::max_value<RealT>(), tol );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(pdf, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1), static_cast<RealT>(1.5) };
    const std::size_t n = sizeof(probs)/sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    // Mathematica: Table[PDF[HyperexponentialDistribution[{0.2, 0.3, 0.5}, {.5, 1.0, 1.5}], x], {x, 0, 4}]
    BOOST_CHECK_CLOSE( boost::math::pdf(dist, static_cast<RealT>(0)), static_cast<RealT>(1.15), tol );
    BOOST_CHECK_CLOSE( boost::math::pdf(dist, static_cast<RealT>(1)), static_cast<RealT>(0.3383645184340184), tol );
    BOOST_CHECK_CLOSE( boost::math::pdf(dist, static_cast<RealT>(2)), static_cast<RealT>(0.11472883036402601), tol );
    BOOST_CHECK_CLOSE( boost::math::pdf(dist, static_cast<RealT>(3)), static_cast<RealT>(0.04558088392888389), tol );
    BOOST_CHECK_CLOSE( boost::math::pdf(dist, static_cast<RealT>(4)), static_cast<RealT>(0.02088728412278129), tol );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(cdf, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs)/sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    // Mathematica: Table[CDF[HyperexponentialDistribution[{0.2, 0.3, 0.5}, {.5, 1.0, 1.5}], x], {x, 0, 4}]
    BOOST_CHECK_CLOSE( boost::math::cdf(dist, static_cast<RealT>(0)), static_cast<RealT>(0), tol );
    BOOST_CHECK_CLOSE( boost::math::cdf(dist, static_cast<RealT>(1)), static_cast<RealT>(0.6567649556318257), tol );
    BOOST_CHECK_CLOSE( boost::math::cdf(dist, static_cast<RealT>(2)), static_cast<RealT>(0.8609299926107957), tol );
    BOOST_CHECK_CLOSE( boost::math::cdf(dist, static_cast<RealT>(3)), static_cast<RealT>(0.9348833491908337), tol );
    BOOST_CHECK_CLOSE( boost::math::cdf(dist, static_cast<RealT>(4)), static_cast<RealT>(0.966198875597724), tol );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(quantile, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs)/sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    // Mathematica: Table[Quantile[HyperexponentialDistribution[{0.2, 0.3, 0.5}, {.5, 1.0, 1.5}], p], {p, {0, 0.6567649556318257, 0.8609299926107957, 0.9348833491908337, 0.966198875597724}}]
    BOOST_CHECK_CLOSE( boost::math::quantile(dist, static_cast<RealT>(0)), static_cast<RealT>(0), tol );
    BOOST_CHECK_CLOSE( boost::math::quantile(dist, static_cast<RealT>(0.6567649556318257)), static_cast<RealT>(1.0000000000000036), tol );
    BOOST_CHECK_CLOSE( boost::math::quantile(dist, static_cast<RealT>(0.8609299926107957)), static_cast<RealT>(1.9999999999999947), tol );
    BOOST_CHECK_CLOSE( boost::math::quantile(dist, static_cast<RealT>(0.9348833491908337)), static_cast<RealT>(3), tol );
    BOOST_CHECK_CLOSE( boost::math::quantile(dist, static_cast<RealT>(0.966198875597724)), static_cast<RealT>(3.9999999999999964), tol );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ccdf, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs)/sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    // Mathematica: Table[SurvivalFunction[HyperexponentialDistribution[{0.2, 0.3, 0.5}, {.5, 1.0, 1.5}], x], {x, 0, 4}]
    BOOST_CHECK_CLOSE( boost::math::cdf(boost::math::complement(dist, static_cast<RealT>(0))), static_cast<RealT>(1), tol );
    BOOST_CHECK_CLOSE( boost::math::cdf(boost::math::complement(dist, static_cast<RealT>(1))), static_cast<RealT>(0.3432350443681743), tol );
    BOOST_CHECK_CLOSE( boost::math::cdf(boost::math::complement(dist, static_cast<RealT>(2))), static_cast<RealT>(0.13907000738920425), tol );
    BOOST_CHECK_CLOSE( boost::math::cdf(boost::math::complement(dist, static_cast<RealT>(3))), static_cast<RealT>(0.0651166508091663), tol );
    BOOST_CHECK_CLOSE( boost::math::cdf(boost::math::complement(dist, static_cast<RealT>(4))), static_cast<RealT>(0.03380112440227598), tol );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(cquantile, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs) / sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    // Mathematica: Table[SurvivalFunction[HyperexponentialDistribution[{0.2, 0.3, 0.5}, {.5, 1.0, 1.5}], p], {p, {1., 0.3432350443681743, 0.13907000738920425, 0.0651166508091663, 0.03380112440227598}}]
    BOOST_CHECK_CLOSE( boost::math::quantile(boost::math::complement(dist, static_cast<RealT>(1))), static_cast<RealT>(0), tol );
    BOOST_CHECK_CLOSE( boost::math::quantile(boost::math::complement(dist, static_cast<RealT>(0.3432350443681743))), static_cast<RealT>(1.0000000000000036), tol );
    BOOST_CHECK_CLOSE( boost::math::quantile(boost::math::complement(dist, static_cast<RealT>(0.13907000738920425))), static_cast<RealT>(1.9999999999999947), tol );
    BOOST_CHECK_CLOSE( boost::math::quantile(boost::math::complement(dist, static_cast<RealT>(0.0651166508091663))), static_cast<RealT>(3), tol );
    BOOST_CHECK_CLOSE( boost::math::quantile(boost::math::complement(dist, static_cast<RealT>(0.03380112440227598))), static_cast<RealT>(3.9999999999999964), tol );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(mean, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs) / sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    // Mathematica: Mean[HyperexponentialDistribution[{0.2, 0.3, 0.5}, {.5, 1.0, 1.5}]]
    BOOST_CHECK_CLOSE( boost::math::mean(dist), static_cast<RealT>(1.0333333333333332), tol );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(variance, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs) / sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    // Mathematica: Mean[HyperexponentialDistribution[{0.2, 0.3, 0.5}, {.5, 1.0, 1.5}]]
    BOOST_CHECK_CLOSE( boost::math::variance(dist), static_cast<RealT>(1.5766666666666673), tol );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(kurtosis, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs) / sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    // Mathematica: Kurtosis[HyperexponentialDistribution[{0.2, 0.3, 0.5}, {.5, 1.0, 1.5}]]
    BOOST_CHECK_CLOSE( boost::math::kurtosis(dist), static_cast<RealT>(19.75073861680871), tol );
    BOOST_CHECK_CLOSE( boost::math::kurtosis_excess(dist), static_cast<RealT>(19.75073861680871)-static_cast<RealT>(3), tol );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(skewness, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs) / sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    // Mathematica: Skewness[HyperexponentialDistribution[{0.2, 0.3, 0.5}, {.5, 1.0, 1.5}]]
    BOOST_CHECK_CLOSE( boost::math::skewness(dist), static_cast<RealT>(3.181138744996378), tol );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(mode, RealT, test_types)
{
    const RealT tol = make_tolerance<RealT>();

    const RealT probs[] = { static_cast<RealT>(0.2L), static_cast<RealT>(0.3L), static_cast<RealT>(0.5L) };
    const RealT rates[] = { static_cast<RealT>(0.5L), static_cast<RealT>(1.0L), static_cast<RealT>(1.5L) };
    const std::size_t n = sizeof(probs) / sizeof(RealT);

    boost::math::hyperexponential_distribution<RealT> dist(probs, probs+n, rates, rates+n);

    BOOST_CHECK_CLOSE( boost::math::mode(dist), static_cast<RealT>(0), tol );
}

