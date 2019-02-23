//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"

BOOST_AUTO_TEST_SUITE(test_autodiff_4)

struct lround_llround_truncl_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::lround;
    using std::llround;
    //using std::truncl; // truncl not supported by boost::multiprecision types.
    constexpr int m = 3;
    const T cx = 3.25;
    auto x = make_fvar<T,m>(cx);
    long yl = lround(x);
    BOOST_REQUIRE(yl == lround(cx));
    long long yll = llround(x);
    BOOST_REQUIRE(yll == llround(cx));
    //long double yld = truncl(x);
    //BOOST_REQUIRE(yld == truncl(cx));
  }
};

BOOST_AUTO_TEST_CASE(lround_llround_truncl)
{
    boost::fusion::for_each(bin_float_types, lround_llround_truncl_test());
    boost::fusion::for_each(multiprecision_float_types, lround_llround_truncl_test());
}

struct multiprecision_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 4000*std::numeric_limits<T>::epsilon(); // percent
    constexpr int Nw=3;
    constexpr int Nx=2;
    constexpr int Ny=4;
    constexpr int Nz=3;
    const auto w = make_fvar<T,Nw>(11);
    const auto x = make_fvar<T,0,Nx>(12);
    const auto y = make_fvar<T,0,0,Ny>(13);
    const auto z = make_fvar<T,0,0,0,Nz>(14);
    const auto v = mixed_partials_f(w,x,y,z); // auto = autodiff_fvar<T,Nw,Nx,Ny,Nz>
    // Calculated from Mathematica symbolic differentiation.
    const T answer = boost::lexical_cast<T>("1976.3196007477977177798818752904187209081211892187"
        "5499076582535951111845769110560421820940516423255314");
    BOOST_REQUIRE_CLOSE(v.derivative(Nw,Nx,Ny,Nz), answer, eps); // Doesn't work for cpp_dec_float
  }
};

BOOST_AUTO_TEST_CASE(multiprecision)
{
    //multiprecision_test()(boost::multiprecision::cpp_bin_float_50());
    boost::fusion::for_each(bin_float_types, multiprecision_test());
}

BOOST_AUTO_TEST_SUITE_END()
