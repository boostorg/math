// Copyright Matt Borland 2025.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/logistic_sigmoid.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "math_unit_test.hpp"
#include <array>
#include <cfloat>
#include <cfenv>

#pragma STDC FENV_ACCESS ON

template <typename RealType>
void test()
{
    const std::array<RealType, 5> x_values = {
        0,
        1,
        1000,
        0.5,
        0.75
    };
    const std::array<RealType, 5> y_values = {
        RealType{1} / 2,
        RealType{0.73105857863000487925115924182183627436514464016505651927636590791904045307},
        RealType{1},
        RealType{0.62245933120185456463890056574550847875327936530891016305943716265854500},
        RealType{0.6791786991753929731596801157765790212342212482195760219829517436805}
    };

    for (std::size_t i = 0; i < x_values.size(); ++i)
    {
        const RealType test_value {boost::math::logistic_sigmoid(x_values[i])};
        BOOST_MATH_IF_CONSTEXPR (std::is_same<RealType, float>::value || std::is_same<RealType, double>::value)
        {
            CHECK_ULP_CLOSE(test_value, y_values[i], 1);
        }
        else
        {
            CHECK_MOLLIFIED_CLOSE(test_value, y_values[i], 1e-15);
        }

        bool fe {false};
        if (std::fetestexcept(FE_OVERFLOW))
        {
            fe = true;
            std::cerr << "FE_OVERFLOW" << std::endl;
        }
        if (std::fetestexcept(FE_UNDERFLOW))
        {
            fe = true;
            std::cerr << "FE_UNDERFLOW" << std::endl;
        }
        if (std::fetestexcept(FE_DIVBYZERO))
        {
            fe = true;
            std::cerr << "FE_DIVBYZERO" << std::endl;
        }
        if (std::fetestexcept(FE_INVALID))
        {
            fe = true;
            std::cerr << "FE_INVALID" << std::endl;
        }

        CHECK_EQUAL(fe, false);
    }
}

int main()
{
    std::feclearexcept(FE_ALL_EXCEPT);
    test<float>();

    std::feclearexcept(FE_ALL_EXCEPT);
    test<double>();

    std::feclearexcept(FE_ALL_EXCEPT);
    test<long double>();

    std::feclearexcept(FE_ALL_EXCEPT);
    test<boost::multiprecision::cpp_bin_float_quad>();

    return boost::math::test::report_errors();
}
