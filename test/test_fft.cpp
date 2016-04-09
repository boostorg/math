//  (C) Copyright Jeremy Murphy 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/config.hpp>
#define BOOST_TEST_MAIN

#include <boost/math/tools/fft.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <vector>

typedef boost::mpl::list<int, long> signed_integral_test_types;
// , boost::multiprecision::cpp_int

using namespace boost;
using namespace boost::math::tools;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_basic, T, signed_integral_test_types)
{
    std::vector<double> a;
    std::copy(make_counting_iterator(1), make_counting_iterator(5), std::back_inserter(a));
    std::vector< std::complex<double> > result;
    result.resize(4);
    fft(a.begin(), a.end(), result.begin());
}
