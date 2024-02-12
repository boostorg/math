// Copyright Hubert Holin 2001.
// Copyright Christopher Kormanyos 2024
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <iomanip>

#include <boost/mpl/list.hpp>
#include <boost/math/octonion.hpp>
#include <boost/core/lightweight_test.hpp>

// test file for octonion.hpp

namespace
{
    template<typename T>
    ::boost::math::octonion<T> index_i_element(int idx)
    {
        return(
            ::boost::math::octonion<T>(
                        (idx == 0) ?
                            static_cast<T>(1) :
                            static_cast<T>(0),
                        (idx == 1) ?
                            static_cast<T>(1) :
                            static_cast<T>(0),
                        (idx == 2) ?
                            static_cast<T>(1) :
                            static_cast<T>(0),
                        (idx == 3) ?
                            static_cast<T>(1) :
                            static_cast<T>(0),
                        (idx == 4) ?
                            static_cast<T>(1) :
                            static_cast<T>(0),
                        (idx == 5) ?
                            static_cast<T>(1) :
                            static_cast<T>(0),
                        (idx == 6) ?
                            static_cast<T>(1) :
                            static_cast<T>(0),
                        (idx == 7) ?
                            static_cast<T>(1) :
                            static_cast<T>(0)
            ));
    }
}

template<class T>
void multiplication_test()
{
    using ::std::numeric_limits;

    using ::boost::math::abs;

    // Testing multiplication.

    const auto one_by_one = ::boost::math::octonion<T>(1,0,0,0,0,0,0,0) * ::boost::math::octonion<T>(1,0,0,0,0,0,0,0);

    const T delta { abs(one_by_one - static_cast<T>(1)) };

    const auto result_mul_one_is_ok = (delta < numeric_limits<T>::epsilon());

    BOOST_TEST(result_mul_one_is_ok);

    for    (int idx = 1; idx < 8; ++idx)
    {
        ::boost::math::octonion<T> toto = index_i_element<T>(idx);

        const T tabs { abs(toto*toto+static_cast<T>(1)) };

        const auto result_mul_toto_is_ok = (tabs < numeric_limits<T>::epsilon());

        BOOST_TEST(result_mul_toto_is_ok);
    }

    {
        const boost::math::octonion<T> lhs(T(1),T(2),T(3),T(4),T(5),T(6),T(7),T(8));
        const boost::math::octonion<T> rhs(T(8),T(7),T(6),T(5),T(4),T(3),T(2),T(1));

        const boost::math::octonion<T> prod = lhs * rhs;

        const boost::math::octonion<T> ctrl(T(-104), T(14), T(12), T(10), T(152), T(42), T(4), T(74));

        BOOST_TEST(prod == ctrl);
    }
}

void octonion_original_manual_test()
{
    // tests for evaluation by humans

    // using default constructor
    ::boost::math::octonion<float>            o0;

    ::boost::math::octonion<float>            oa[2];

    // using constructor "O seen as R^8"
    ::boost::math::octonion<float>            o1(1,2,3,4,5,6,7,8);

    ::std::complex<double>                    c0(9,10);

    // using constructor "O seen as C^4"
    ::boost::math::octonion<double>            o2(c0);

    ::boost::math::quaternion<long double>    q0(11,12,13,14);

    // using constructor "O seen as H^2"
    ::boost::math::octonion<long double>      o3(q0);

    // using UNtemplated copy constructor
    ::boost::math::octonion<float>            o4(o1);

    // using templated copy constructor
    ::boost::math::octonion<long double>      o5(o2);

    // using UNtemplated assignment operator
    o5 = o3;
    oa[0] = o0;

    // using templated assignment operator
    o5 = o2;
    oa[1] = o5;

    float                                     f0(15);

    // using converting assignment operator
    o0 = f0;

    // using converting assignment operator
    o2 = c0;

    // using converting assignment operator
    o5 = q0;

    // using += (const T &)
    o4 += f0;

    // using == (const octonion<T> &,const octonion<T> &)
    BOOST_TEST(o0 != o4);

    // using += (const ::std::complex<T> &)
    o2 += c0;

    // using == (const ::boost::math::quaternion<T> &, const octonion<T> &)
    BOOST_TEST(q0 == o3);

    // using += (const ::boost::math::quaternion<T> &)
    o3 += q0;

    BOOST_TEST(2 * q0 == o3);

    // using += (const quaternion<X> &)
    o5 += o4;

    // using -= (const T &)
    o1 -= f0;

    // using -= (const ::std::complex<T> &)
    o2 -= c0;

    // using -= (const ::boost::math::quaternion<T> &)
    o5 -= q0;

    // using -= (const octonion<X> &)
    o3 -= o4;

    // using == (const ::std::complex<T> &, const octonion<T> &)
    BOOST_TEST(c0 == o2);

    // using == (const octonion<T> &, const ::std::complex<T> &)
    BOOST_TEST(o2 == c0);

    double                                    d0(16);
    ::std::complex<double>                    c1(17,18);
    ::boost::math::quaternion<double>         q1(19,20,21,22);

    // using *= (const T &)
    o2 *= d0;

    // using *= (const ::std::complex<T> &)
    o2 *= c1;

    // using *= (const ::boost::math::quaternion<T> &)
    o2 *= q1;

    // using *= (const octonion<X> &)
    o2 *= o4;

    long double                               l0(23);
    ::std::complex<long double>               c2(24,25);

    // using /= (const T &)
    o5 /= l0;

    // using /= (const ::std::complex<T> &)
    o5 /= c2;

    // using /= (const quaternion<X> &)
    o5 /= q0;

    // using /= (const octonion<X> &)
    o5 /= o5;

    // using + (const T &, const octonion<T> &)
    ::boost::math::octonion<float>            o6 = f0+o0;

    // using + (const octonion<T> &, const T &)
    ::boost::math::octonion<float>            o7 = o0+f0;

    // using + (const ::std::complex<T> &, const quaternion<T> &)
    ::boost::math::octonion<double>           o8 = c0+o2;

    // using + (const octonion<T> &, const ::std::complex<T> &)
    ::boost::math::octonion<double>           o9 = o2+c0;

    // using + (const ::boost::math::quaternion<T>, const octonion<T> &)
    ::boost::math::octonion<long double>      o10 = q0+o3;

    // using + (const octonion<T> &, const ::boost::math::quaternion<T> &)
    ::boost::math::octonion<long double>      o11 = o3+q0;

    // using + (const quaternion<T> &,const quaternion<T> &)
    ::boost::math::octonion<float>            o12 = o0+o4;

    // using - (const T &, const octonion<T> &)
    o6 = f0-o0;

    // using - (const octonion<T> &, const T &)
    o7 = o0-f0;

    // using - (const ::std::complex<T> &, const octonion<T> &)
    o8 = c0-o2;

    // using - (const octonion<T> &, const ::std::complex<T> &)
    o9 = o2-c0;

    // using - (const quaternion<T> &,const octonion<T> &)
    o10 = q0-o3;

    // using - (const octonion<T> &,const quaternion<T> &)
    o11 = o3-q0;

    // using - (const octonion<T> &,const octonion<T> &)
    o12 = o0-o4;

    // using * (const T &, const octonion<T> &)
    o6 = f0*o0;

    // using * (const octonion<T> &, const T &)
    o7 = o0*f0;

    // using * (const ::std::complex<T> &, const octonion<T> &)
    o8 = c0*o2;

    // using * (const octonion<T> &, const ::std::complex<T> &)
    o9 = o2*c0;

    // using * (const quaternion<T> &,const octonion<T> &)
    o10 = q0*o3;

    // using * (const octonion<T> &,const quaternion<T> &)
    o11 = o3*q0;

    // using * (const octonion<T> &,const octonion<T> &)
    o12 = o0*o4;

    // using / (const T &, const octonion<T> &)
    o6 = f0/o0;

    // using / (const octonion<T> &, const T &)
    o7 = o0/f0;

    // using / (const ::std::complex<T> &, const octonion<T> &)
    o8 = c0/o2;

    // using / (const octonion<T> &, const ::std::complex<T> &)
    o9 = o2/c0;

    // using / (const ::boost::math::quaternion<T> &, const octonion<T> &)
    o10 = q0/o3;

    // using / (const octonion<T> &, const ::boost::math::quaternion<T> &)
    o11 = o3/q0;

    // using / (const octonion<T> &,const octonion<T> &)
    o12 = o0/o4;

    // using + (const octonion<T> &)
    o4 = +o0;

    // using == (const T &, const octonion<T> &)
    BOOST_TEST(f0 == o0);

    // using == (const octonion<T> &, const T &)
    BOOST_TEST(o0 == f0);

    // using - (const octonion<T> &)
    o0 = -o4;

    // using != (const T &, const octonion<T> &)
    BOOST_TEST(f0 != o0);

    // using != (const octonion<T> &, const T &)
    BOOST_TEST(o0 != f0);

    // using != (const ::std::complex<T> &, const octonion<T> &)
    BOOST_TEST(c0 != o2);

    // using != (const octonion<T> &, const ::std::complex<T> &)
    BOOST_TEST(o2 != c0);

    // using != (const ::boost::math::quaternion<T> &, const octonion<T> &)
    BOOST_TEST(q0 != o3);

    // using != (const octonion<T> &, const ::boost::math::quaternion<T> &)
    BOOST_TEST(o3 != q0);

    // using != (const octonion<T> &,const octonion<T> &)
    BOOST_TEST(o0 != o4);
}

template <class T>
void exp_test()
{
    using ::std::numeric_limits;

    using ::std::atan;

    using ::boost::math::abs;

    // Testing exp.

    for(int idx = 1; idx < 8; ++idx)
    {
        ::boost::math::octonion<T> toto =
            static_cast<T>(4)*atan(static_cast<T>(1))*index_i_element<T>(idx);

        const T tabs { abs(exp(toto)+static_cast<T>(1)) };

        const auto result_exp_toto_is_ok = (tabs < 2*numeric_limits<T>::epsilon());

        BOOST_TEST(result_exp_toto_is_ok);
    }
}

auto main() -> int
{
  multiplication_test<float>();
  multiplication_test<double>();

  exp_test<float>();
  exp_test<double>();

  octonion_original_manual_test();

  const auto result_is_ok = (boost::report_errors() == 0);

  return (result_is_ok ? 0 : -1);
}
