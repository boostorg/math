// Copyright Matt Borland 2025.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/logit.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "math_unit_test.hpp"
#include <array>
#include <limits>
#include <cfenv>

#pragma STDC FENV_ACCESS ON

template <typename RealType>
void test()
{
    const std::array<RealType, 5> x_values = {
        0.01,
        0.24,
        0.5,
        0.75,
        0.995,
    };
    const std::array<RealType, 5> y_values = {
        RealType{-4.595119850134589926852434051810180709116687969582916078687956376405},
        RealType{-1.15267950993838545919655007350715126451438856911612411268258589327840479},
        RealType{0},
        RealType{1.09861228866810969139524523692252570464749055782274945173469433363749429},
        RealType{5.2933048247244923954101212918685372018911052805694724989064609879440992}
    };

    for (std::size_t i = 0; i < x_values.size(); ++i)
    {
        const RealType test_value {boost::math::logit(x_values[i])};

        BOOST_MATH_IF_CONSTEXPR (std::is_same<RealType, float>::value || std::is_same<RealType, double>::value)
        {
            CHECK_ULP_CLOSE(test_value, y_values[i], 5);
        }
        else
        {
            CHECK_MOLLIFIED_CLOSE(test_value, y_values[i], 1e-15);
        }

        bool fe {false};
        if (std::fetestexcept(FE_OVERFLOW))
        {
            fe = true;                                  // LCOV_EXCL_LINE
            std::cerr << "FE_OVERFLOW" << std::endl;    // LCOV_EXCL_LINE
        }
        if (std::fetestexcept(FE_UNDERFLOW))
        {
            fe = true;                                  // LCOV_EXCL_LINE
            std::cerr << "FE_UNDERFLOW" << std::endl;   // LCOV_EXCL_LINE
        }
        if (std::fetestexcept(FE_DIVBYZERO))
        {
            fe = true;                                  // LCOV_EXCL_LINE
            std::cerr << "FE_DIVBYZERO" << std::endl;   // LCOV_EXCL_LINE
        }
        if (std::fetestexcept(FE_INVALID))
        {
            fe = true;                                  // LCOV_EXCL_LINE
            std::cerr << "FE_INVALID" << std::endl;     // LCOV_EXCL_LINE
        }

        CHECK_EQUAL(fe, false);
    }

    #if defined(_CPPUNWIND) || defined(__EXCEPTIONS)

    BOOST_MATH_IF_CONSTEXPR (std::is_arithmetic<RealType>::value)
    {
        bool thrown {false};
        try
        {
            boost::math::logit(std::numeric_limits<RealType>::denorm_min());
        }
        catch (...)
        {
            thrown = true;
        }

        CHECK_EQUAL(thrown, true);
    }

    #endif // Exceptional environments
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
