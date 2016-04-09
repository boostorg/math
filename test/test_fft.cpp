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

typedef boost::mpl::list<int> signed_integral_test_types;
// , boost::multiprecision::cpp_int

using namespace boost;
using namespace boost::math::tools;

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


BOOST_AUTO_TEST_CASE_TEMPLATE(test_bit_reversed_permute, T, signed_integral_test_types)
{
    std::vector<int> a;
    std::copy(make_counting_iterator(0), make_counting_iterator(16), std::back_inserter(a));
    bit_reversed_permute(a.begin(), a.end());
    boost::array<int, 16> b = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};
    BOOST_CHECK_EQUAL_COLLECTIONS(boost::begin(a), boost::end(a), boost::begin(b), boost::end(b));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_iterative_fft, T, signed_integral_test_types)
{
    std::vector<double> a;
    std::copy(make_counting_iterator(0), make_counting_iterator(16), std::back_inserter(a));
    std::vector< std::complex<double> > result;
    result.resize(16);
    iterative_fft(boost::begin(a), boost::end(a), boost::begin(result));
    
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
    iterative_fft(boost::begin(a), boost::end(a), boost::begin(iterative_result));
    recursive_fft(a.begin(), a.end(), recursive_result.begin());
    BOOST_CHECK_EQUAL_COLLECTIONS(iterative_result.begin(), iterative_result.end(), recursive_result.begin(), recursive_result.end());
}
