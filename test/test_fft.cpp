//  (C) Copyright Jeremy Murphy 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/config.hpp>
#define BOOST_TEST_MAIN

#include <boost/array.hpp>
#include <boost/math/tools/fft.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/algorithm.hpp>

#include <gsl/gsl_fft_real.h>

#include <vector>
#include <iostream>
#include <iterator>
#include <functional>

typedef boost::mpl::list<int> signed_integral_test_types;
// , boost::multiprecision::cpp_int

using namespace boost;
using namespace boost::math::tools;
using boost::math::tools::detail::identity;

typedef std::complex<double> complex_t;

// TODO: Fixtures.

BOOST_AUTO_TEST_CASE_TEMPLATE(test_recursive_fft, T, signed_integral_test_types)
{
    std::vector<double> a;
    std::copy(make_counting_iterator(0), make_counting_iterator(16), std::back_inserter(a));
    std::vector< std::complex<double> > result;
    result.resize(16);
    recursive_fft(a.begin(), a.end(), result.begin());
    // TODO: Compare against benchmark result (GSL).
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_iterative_fft_radix2, T, signed_integral_test_types)
{
    std::vector<double> a;
    std::copy(make_counting_iterator(0), make_counting_iterator(16), std::back_inserter(a));
    std::vector< std::complex<double> > result;
    result.resize(16);
    iterative_fft_radix2(boost::begin(a), boost::end(a), boost::begin(result), identity<complex_t>());
    
    gsl_fft_real_radix2_transform(a.data(), 1, a.size());
   
    // TODO: Compare against benchmark result (GSL).
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_fft_agreement, T, signed_integral_test_types)
{
    // Tests whether the two implementations produce the same result.
    std::vector<double> a;
    std::copy(make_counting_iterator(0), make_counting_iterator(16), std::back_inserter(a));
    std::vector< std::complex<double> > recursive_result;
    std::vector< std::complex<double> > iterative_result;
    recursive_result.resize(16);
    iterative_result.resize(16);
    iterative_fft_radix2(boost::begin(a), boost::end(a), boost::begin(iterative_result), identity<complex_t>());
    recursive_fft(a.begin(), a.end(), recursive_result.begin());
    BOOST_CHECK_EQUAL_COLLECTIONS(iterative_result.begin(), iterative_result.end(), recursive_result.begin(), recursive_result.end());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_fft_involution, T, signed_integral_test_types)
{
    std::vector<double> a;
    std::copy(make_counting_iterator(0), make_counting_iterator(16), std::back_inserter(a));
    std::vector< std::complex<double> > tmp;
    tmp.resize(16);
    iterative_fft_radix2_forward(a.begin(), a.end(), tmp.begin());
    std::vector<complex_t> b;
    b.resize(16);
    iterative_fft_radix2_inverse(tmp.begin(), tmp.end(), b.begin());
    for (std::size_t i = 0; i < a.size(); i++)
        BOOST_CHECK_EQUAL(a[0], b[0].real());
}
