//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"

BOOST_AUTO_TEST_SUITE(test_autodiff_4)

struct round_and_trunc_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::round;
    using std::trunc;
    constexpr int m = 3;
    constexpr float cx = 3.25;
    auto x = make_fvar<T,m>(cx);
    auto y = round(x);
    BOOST_REQUIRE(y.derivative(0) == round(cx));
    BOOST_REQUIRE(y.derivative(1) == 0.0);
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
    y = trunc(x);
    BOOST_REQUIRE(y.derivative(0) == trunc(cx));
    BOOST_REQUIRE(y.derivative(1) == 0.0);
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(round_and_trunc)
{
    boost::fusion::for_each(bin_float_types, round_and_trunc_test());
    boost::fusion::for_each(multiprecision_float_types, round_and_trunc_test());
}

struct iround_and_itrunc_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using namespace boost::math;
    constexpr int m = 3;
    constexpr float cx = 3.25;
    auto x = make_fvar<T,m>(cx);
    int y = iround(x);
    BOOST_REQUIRE(y == iround(cx));
    y = itrunc(x);
    BOOST_REQUIRE(y == itrunc(cx));
  }
};

BOOST_AUTO_TEST_CASE(iround_and_itrunc)
{
    boost::fusion::for_each(bin_float_types, iround_and_itrunc_test());
    boost::fusion::for_each(multiprecision_float_types, iround_and_itrunc_test());
}

struct lambert_w0_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 1000*std::numeric_limits<T>::epsilon(); // percent
    constexpr int m = 10;
    const T cx = 3;
    // Mathematica: N[Table[D[ProductLog[x], {x, n}], {n, 0, 10}] /. x -> 3, 52]
    const char* const answers[m+1] {
        "1.049908894964039959988697070552897904589466943706341",
        "0.1707244807388472968312949774415522047470762509741737",
        "-0.04336545501146252734105411312976167858858970875797718",
        "0.02321456264324789334313200360870492961288748451791104",
        "-0.01909049778427783072663170526188353869136655225133878",
        "0.02122935002563637629500975949987796094687564718834156",
        "-0.02979093848448877259041971538394953658978044986784643",
        "0.05051290266216717699803334605370337985567016837482099",
        "-0.1004503154972645060971099914384090562800544486549660",
        "0.2292464437392250211967939182075930820454464472006425",
        "-0.5905839053125614593682763387470654123192290838719517"};
    auto x = make_fvar<T,m>(cx);
    auto y = lambert_w0(x);
    for (int i=0 ; i<=m ; ++i)
    {
        const T answer = boost::lexical_cast<T>(answers[i]);
        BOOST_REQUIRE_CLOSE(y.derivative(i), answer, eps);
    }
    //const T cx0 = -1 / boost::math::constants::e<T>();
    //auto edge = lambert_w0(make_fvar<T,m>(cx0));
    //std::cout << "edge = " << edge << std::endl;
    //edge = depth(1)(-1,inf,-inf,inf,-inf,inf,-inf,inf,-inf,inf,-inf)
    //edge = depth(1)(-1,inf,-inf,inf,-inf,inf,-inf,inf,-inf,inf,-inf)
    //edge = depth(1)(-1,3.68935e+19,-9.23687e+57,4.62519e+96,-2.89497e+135,2.02945e+174,-1.52431e+213,1.19943e+252,-9.75959e+290,8.14489e+329,-6.93329e+368)
  }
};

BOOST_AUTO_TEST_CASE(lambert_w0_test)
{
    boost::fusion::for_each(bin_float_types, lambert_w0_test_test());
    boost::fusion::for_each(multiprecision_float_types, lambert_w0_test_test());
}

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
    const T eps = 600*std::numeric_limits<T>::epsilon(); // percent
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
    // BOOST_REQUIRE_CLOSE(v.derivative(Nw,Nx,Ny,Nz), answer, eps); // Doesn't work for cpp_dec_float
    using std::fabs;
    const double relative_error = static_cast<double>(fabs(v.derivative(Nw,Nx,Ny,Nz)/answer-1));
    BOOST_REQUIRE(relative_error < eps);
  }
};

BOOST_AUTO_TEST_CASE(multiprecision)
{
    //multiprecision_test()(boost::multiprecision::cpp_bin_float_50());
    boost::fusion::for_each(bin_float_types, multiprecision_test());
}

BOOST_AUTO_TEST_SUITE_END()
