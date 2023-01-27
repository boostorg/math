// Copyright Matt Borland, 2023
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <climits>
#include <cmath>
#include <stdexcept>
#include <boost/math/special_functions/powm1.hpp>
#include "math_unit_test.hpp"

void test(double x, double y, double p)
{
    double z;
    try {
        z = boost::math::powm1(x, y);
    } catch (const std::domain_error& e) {
        z = NAN;
    } catch (const std::overflow_error& e) {
        
        // See: https://en.cppreference.com/w/cpp/numeric/math/pow
        if (x > 0) {
            if (y < 0) {
                z = 0;
            }
            else if (y == 0) {
                z = 1;
            }
            else {
                z = INFINITY;
            }
        }
        else if (x == 0)
        {
            z = INFINITY;
        }
        else {
            if (y < 0) {
                if (std::fmod(y, 2) == 0) {
                    z = 0;
                }
                else {
                    z = -0;
                }
            }
            else if (y == 0) {
                z = 1;
            }
            else {
                if (std::fmod(y, 2) == 0) {
                    z = INFINITY;
                }
                else {
                    z = -INFINITY;
                }
            }
        }
    } catch (const std::underflow_error& e) {
        z = 0;
    } catch (...) {
        z = NAN;
    }
    
    CHECK_EQUAL(p, z);
}

int main(void)
{
    test(0.0, 0.0, 0.0);
    test(0.0, -1.5, INFINITY);
    test(0.0, 1.75, -1.0);
    test(-1.5, 2.0, 1.25);
    test(-1.5, 3.0, -4.375);
    test(NAN, 0.0, 0.0);
    test(1.0, NAN, 0.0);
    test(1.0, INFINITY, 0.0);
    test(1.0, -INFINITY, 0.0);
    test(INFINITY, 7.5, INFINITY);
    test(INFINITY, -7.5, -1.0);
    test(3.25, INFINITY, INFINITY);
    test(INFINITY, INFINITY, INFINITY);
    test(INFINITY, -INFINITY, -1.0);
    test(INFINITY, 0.0, 0.0);
    test(-INFINITY, 0.0, 0.0);
    test(-INFINITY, 2.0, INFINITY);
    test(-INFINITY, 3.0, -INFINITY);
    test(-1.0, std::pow(2, 53) - 1, -2.0);

    return boost::math::test::report_errors();
}
