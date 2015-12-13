// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Copyright Jeremy W. Murphy 2015.

// This file is written to be included from a Quickbook .qbk document.
// It can be compiled by the C++ compiler, and run. Any output can
// also be added here as comment or included or pasted in elsewhere.
// Caution: this file contains Quickbook markup as well as code
// and comments: don't change any of the special comment markups!


//[polynomial_arithmetic1
/*
Include the essential polynomial header and others to make the example.
*/
#include <boost/math/tools/polynomial.hpp>
#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assert.hpp>

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <string>
#include <utility>

using std::string;
using std::exception;
using std::cout;
using std::abs;
using std::pair;

using namespace boost::math;
using namespace boost::math::tools;
using boost::lexical_cast;

//] [/polynomial_arithmetic1]


template <typename T>
string sign_str(T const &x)
{
    return x < 0 ? "-" : "+";
}

template <typename T>
string inner_coefficient(T const &x)
{
    string result(" " + sign_str(x) + " ");
    if (abs(x) != T(1))
        result += lexical_cast<string>(abs(x));
    return result;
}

template <typename T>
string pretty_format(polynomial<T> const &a)
{
    string result;
    if (a.size() == 0)
        result += lexical_cast<string>(T(0));
    else
    {
        // First one is a special case as it may need unary negate.
        unsigned i = a.size() - 1;
        if (a[i] < 0)
            result += "-";
        if (abs(a[i]) != T(1))
            result += lexical_cast<string>(abs(a[i]));
        
        if (i > 0)
        {
            result += "x";
            if (i > 1)
            {
                result += "^" + lexical_cast<string>(i);
                i--;
                for (; i != 1; i--)
                    result += inner_coefficient(a[i]) + "x^" + lexical_cast<string>(i);
                
                result += inner_coefficient(a[i]) + "x";
            }
            i--;
            
            result += " " + sign_str(a[i]) + " " + lexical_cast<string>(abs(a[i]));
        }
    }
    return result;
}


int main()
{
    cout << "Example: Polynomial arithmetic.\n\n";
    
    try
    {
//[polynomial_arithmetic2
/*
Store the coefficients in a convenient way to access them, then create some
polynomials using construction from an iterator range.

Note that although we might conventionally write a polynomial from left to right 
in descending order of degree, they are stored in ascending order of degree.

read/write for humans:    3x^3 - 4x^2 - 6x + 10
Boost polynomial storage: [ 10, -6, -4, 3 ]
*/
        boost::array<double, 4> const d3a = {{10, -6, -4, 3}};
        boost::array<double, 2> const d1a = {{-2, 1}};
        polynomial<double> const a(d3a.begin(), d3a.end());
        polynomial<double> const b(d1a.begin(), d1a.end());
// pretty_format() is doing the conversion from Boost storage to human notation.
        cout << "a = " << pretty_format(a) << "\nb = " << pretty_format(b) << "\n\n";
// Now we can do arithmetic with the usual infix operators.
        polynomial<double> s = a + b;
        cout << "a + b = " << pretty_format(s) << "\n";
        polynomial<double> d = a - b;
        cout << "a - b = " << pretty_format(d) << "\n";
        polynomial<double> p = a * b;
        cout << "a * b = " << pretty_format(p) << "\n";
        polynomial<double> q = a / b;
        cout << "a / b = " << pretty_format(q) << "\n";
        polynomial<double> r = a % b;
        cout << "a % b = " << pretty_format(r) << "\n";
//] [/polynomial_arithmetic2]
        
//[polynomial_arithmetic3
/*
Division is a special case where you can calculate two for the price of one.
Actually, quotient and remainder are always calculated together due to the nature
of the algorithm: the infix operators return one result and throw the other
away.

If you are doing a lot of division and want both the quotient and remainder then
you don't want to do twice the work necessary. In that case you can call the
underlying function, quotient_remainder, to get both results together.
*/
        pair< polynomial<double>, polynomial<double> > result;
        result = quotient_remainder(a, b);
// Reassure ourselves that the result is the same.
        BOOST_ASSERT(result.first == q);
        BOOST_ASSERT(result.second == r);
    }
    catch (exception const &e)
    {
        cout << "\nMessage from thrown exception was:\n   " << e.what() << "\n";
    }
//] [/polynomial_arithmetic3]
}
