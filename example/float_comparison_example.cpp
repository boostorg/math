//!file
//! \brief floating-point comparison from Boost.Test
// Copyright Paul A. Bristow 2015.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Note that this file contains Quickbook mark-up as well as code
// and comments, don't change any of the special comment mark-ups!

#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/next.hpp>

#include <iostream>
#include <limits> // for std::numeric_limits<T>::epsilon().

int main()
{
  std::cout << "Compare floats using Boost.Test functions/classes" << std::endl;

//[compare_floats_synopsis
/*`

  namespace boost { namespace math {
  namespace fpc { // Note floating-point comparison namespace.

  template<typename FPT1, typename FPT2, typename ToleranceType>
  bool is_close_to( FPT1 left, FPT2 right, ToleranceType tolerance );
  // Test if two values are close enough,
  // (using the default FPC_STRONG or 'essentially equal' criterion).

  enum strength
  {
    FPC_STRONG, // "Very close" "essentially equal" - Knuth equation 1' in docs (default).
    FPC_WEAK    // "Close enough" "approximately equal" - equation 2' in docs.
  };

  template<typename ToleranceType>
  explicit close_at_tolerance(ToleranceType tolerance, fpc::strength fpc_strength = FPC_STRONG );

Comparisons are most simply made using the function `is_close_to`.

There is also a templated class `close_at_tolerance` that can be convenient
for multiple tests with the same tolerance and strength.

(These are used by the popular MACRO versions in Boost.Test like `BOOST_CHECK_CLOSE`).

For most applications, the default strength parameter can be left at the default 'strong'.

The `Tolerance_type` is the same as floating-point type `FPT`, often a built-in type
like `float`, `double` or `long double`,
but also __multiprecision types like __cpp_bin_float or __cpp_dec_float.

The constructor sets the [*fractional] tolerance and the equality strength.

Two member functions allow access to the chosen tolerance and strength.

  FPT fraction_tolerance() const;
  strength strength() const; // weak or strong.

the `operator()` functor carries out the comparison,
and returns `true` if ['essentially equal] else `false`.

  bool operator()(FPT left, FPT right) const; // true if close or 'equal'.

Comparison tolerances can be very small, near the
[@http://en.wikipedia.org/wiki/Machine_epsilon machine epsilon]
or [@https://en.wikipedia.org/wiki/Unit_in_the_last_place Unit in Last Place (ULP)],
typically for measuring 'computational' noise from multiple rounding or iteration,
or can be a much bigger value like 0.01 (equivalent to a 1% tolerance),
typically from measurement uncertainty.

After ([*but not before]) a comparison of values ['u] and ['v]
has been made by a call of the functor `operator()`, the access function

  FPT failed_fraction() const;

returns the fraction

__spaces ['abs(u-v) / abs(v)]  or ['abs(u-v) / abs(u)]

that failed the test.

*/
//] [/compare_floats_synopsis]

//[compare_floats_using
/*`Some using statements will ensure that the classes, functions and enums are accessible.
*/

  using namespace boost::math::fpc;

//`or

  using boost::math::fpc::close_at_tolerance;
  using boost::math::fpc::small_with_tolerance;
  using boost::math::fpc::is_close_to;
  using boost::math::fpc::is_small;
  using boost::math::fpc::FPC_STRONG;
  using boost::math::fpc::FPC_WEAK;

//] [/compare_floats_using]


//[compare_floats_example_1
/*`The following examples display values with all possibly significant digits.
Newer compilers should provide `std::numeric_limitsFPT>::max_digits10`
for this purpose, and here we use `float` precision where `max_digits10` = 9
to avoid displaying a distracting number of decimal digits.

[note Older compilers can use this formula to calculate `max_digits10`
from `std::numeric_limits<FPT>::digits10`:[br]
__spaces `int max_digits10 = 2 + std::numeric_limits<FPT>::digits10 * 3010/10000;`
] [/note]

One can set the display including all trailing zeros
(helpful for this example to show all potentially significant digits),
and also to display `bool` values as words rather than integers:
*/
  std::cout.precision(std::numeric_limits<float>::max_digits10);
  std::cout << std::boolalpha << std::showpoint << std::endl;

//] [/compare_floats_example_1]

//[compare_floats_example_2]
/*`
When comparing values that are ['quite close] or ['approximately equal],
it is convenient to use the appropriate `epsilon` for the floating-point type `FPT`, here, for example, `float`:
*/

  float epsilon = std::numeric_limits<float>::epsilon();
  std::cout << "float epsilon = " << epsilon << std::endl; // +1.1920929e-007

//] [/compare_floats_example_2]

//[compare_floats_example_3
//`The simplest use is to compare two values with a tolerance thus:

  bool is_close = is_close_to(1.F, 1.F + epsilon, epsilon); // One epsilon apart is close enough.
  std::cout << "is_close_to(1.F, 1.F + epsilon, epsilon); is " << is_close << std::endl; // true

  is_close = is_close_to(1.F, 1.F + 2 * epsilon, epsilon); // Two epsilon apart isn't close enough.
  std::cout << "is_close_to(1.F, 1.F + epsilon, epsilon); is " << is_close << std::endl; // false

/*`
[note The type FPT of the tolerance and the type of the values [*must match].

So `is_close(0.1F, 1., 1.)` will fail to compile because "template parameter 'FPT' is ambiguous".
Always provide the same type, using `static_cast<FPT>` if necessary.]
*/


/*`An instance of class `close_at_tolerance` is more convenient
when multiple tests with the same conditions are planned.
A class that stores a tolerance of three epsilon (and the default ['strong] test) is:
*/

  close_at_tolerance<float> three_rounds(3 * epsilon); // 'strong' by default.

//`and we can confirm these settings:

  std::cout << "fraction_tolerance = "
    << three_rounds.fraction_tolerance()
    << std::endl; // +3.57627869e-007
  std::cout << "strength = "
    << (three_rounds.strength() == FPC_STRONG ? "strong" : "weak")
    << std::endl; // strong

//`To start, let us use two values that are truly equal (having identical bit patterns)

  float a = 1.23456789F;
  float b = 1.23456789F;

//`and make a comparison using our 3*epsilon `three_rounds` functor:

  bool close = three_rounds(a, b);
  std::cout << "three_rounds(a, b) = " << close << std::endl; // true

//`Unsurprisingly, the result is true, and the failed fraction is zero.

  std::cout << "failed_fraction = " << three_rounds.failed_fraction() << std::endl;

/*`To get some nearby values, it is convenient to use the Boost.Math __next_float functions,
for which we need an include

  #include <boost/math/special_functions/next.hpp>

and some using declarations:
*/

  using boost::math::float_next;
  using boost::math::float_prior;
  using boost::math::nextafter;
  using boost::math::float_distance;

//`To add a few __ulp to one value:
  b = float_next(a); // Add just one ULP to a.
  b = float_next(b); // Add another one ULP.
  b = float_next(b); // Add another one ULP.
  // 3 epsilon would pass.
  b = float_next(b); // Add another one ULP.

//`and repeat our comparison:

  close = three_rounds(a, b);
  std::cout << "three_rounds(a, b) = " << close << std::endl; // false
  std::cout << "failed_fraction = " << three_rounds.failed_fraction()
    << std::endl;  // abs(u-v) / abs(v) = 3.86237957e-007

//`We can also 'measure' the number of bits different using the `float_distance` function:

  std::cout << "float_distance = " << float_distance(a, b) << std::endl; // 4

/*`Now consider two values that are much further apart
than one might expect from ['computational noise],
perhaps the result of two measurements of some physical property like length
where an uncertainty of a percent or so might be expected.
*/
  float fp1 = 0.01000F;
  float fp2 = 0.01001F; // Slightly different.

  float tolerance = 0.0001F;

  close_at_tolerance<float> strong(epsilon); // Default is strong.
  bool rs = strong(fp1, fp2);
  std::cout << "strong(fp1, fp2) is " << rs << std::endl;

//`Or we could contrast using the ['weak] criterion:
  close_at_tolerance<float> weak(epsilon, FPC_WEAK); // Explicitly weak.
  bool rw = weak(fp1, fp2); //
  std::cout << "weak(fp1, fp2) is " << rw << std::endl;

//`We can also construct, setting tolerance and strength, and compare in one statement:

  std::cout << a << " #= " << b << " is "
    << close_at_tolerance<float>(epsilon, FPC_STRONG)(a, b) << std::endl;
  std::cout << a << " ~= " << b << " is "
    << close_at_tolerance<float>(epsilon, FPC_WEAK)(a, b) << std::endl;

//`but this has little advantage over using function `is_close_to` directly.

//] [/compare_floats_example_3]


/*When the floating-point values become very small and near zero, using
//a relative test becomes unhelpful because one is dividing by zero or tiny,

//Instead, an absolute test is needed, comparing one (or usually both) values with zero,
//using a tolerance.
//This is provided by the `small_with_tolerance` class and `is_small` function.

  namespace boost {
  namespace math {
  namespace fpc {


  template<typename FPT>
  class small_with_tolerance
  {
  public:
  // Public typedefs.
  typedef bool result_type;

  // Constructor.
  explicit small_with_tolerance(FPT tolerance); // tolerance >= 0

  // Functor
  bool operator()(FPT value) const; // return true if <= absolute tolerance (near zero).
  };

  template<typename FPT>
  bool
  is_small(FPT value, FPT tolerance); // return true if value <= absolute tolerance (near zero).

  }}} // namespaces.

/*`
[note The type FPT of the tolerance and the type of the value [*must match].

So `is_small(0.1F, 0.000001)` will fail to compile because "template parameter 'FPT' is ambiguous".
Always provide the same type, using `static_cast<FPT>` if necessary.]

A few values near zero are tested with varying tolerance below.
*/
//[compare_floats_small_1

  float c = 0;
  std::cout << "0 is_small " << is_small(c, epsilon) << std::endl; // true

  c = std::numeric_limits<float>::denorm_min(); // 1.40129846e-045
  std::cout << "denorm_ min =" << c << ", is_small is " << is_small(c, epsilon) << std::endl; // true

  c = std::numeric_limits<float>::min(); // 1.17549435e-038
  std::cout << "min = " << c << ", is_small is " << is_small(c, epsilon) << std::endl; // true

  c = 1 * epsilon; // 1.19209290e-007
  std::cout << "epsilon = " << c << ", is_small is " << is_small(c, epsilon) << std::endl; // false

  c = 1 * epsilon; // 1.19209290e-007
  std::cout << "2 epsilon = " << c << ", is_small is " << is_small(c, 2 * epsilon) << std::endl; // true

  c = 2 * epsilon; //2.38418579e-007
  std::cout << "4 epsilon = " << c << ", is_small is " << is_small(c, 2 * epsilon) << std::endl; // false

  c = 0.00001F;
  std::cout << "0.00001 = " << c << ", is_small is " << is_small(c, 0.0001F) << std::endl; // true

  c = -0.00001F;
  std::cout << "0.00001 = " << c << ", is_small is " << is_small(c, 0.0001F) << std::endl; // true

/*`Using the class `small_with_tolerance` allows storage of the tolerance,
convenient if you make repeated tests with the same tolerance.
*/

  small_with_tolerance<float>my_test(0.01F);

  std::cout << "my_test(0.001F) is " << my_test(0.001F) << std::endl; // true
  std::cout << "my_test(0.001F) is " << my_test(0.01F) << std::endl; // false

  //] [/compare_floats_small_1]

  return 0;
}  // int main()

/*

Example output is:

//[compare_floats_output
Compare floats using Boost.Test functions/classes

float epsilon = 1.19209290e-007
is_close_to(1.F, 1.F + epsilon, epsilon); is true
is_close_to(1.F, 1.F + epsilon, epsilon); is false
fraction_tolerance = 3.57627869e-007
strength = strong
three_rounds(a, b) = true
failed_fraction = 0.000000000
three_rounds(a, b) = false
failed_fraction = 3.86237957e-007
float_distance = 4.00000000
strong(fp1, fp2) is false
weak(fp1, fp2) is false
1.23456788 #= 1.23456836 is false
1.23456788 ~= 1.23456836 is false
0 is_small true
denorm_ min =1.40129846e-045, is_small is true
min = 1.17549435e-038, is_small is true
epsilon = 1.19209290e-007, is_small is false
2 epsilon = 1.19209290e-007, is_small is true
4 epsilon = 2.38418579e-007, is_small is false
0.00001 = 9.99999975e-006, is_small is true
0.00001 = -9.99999975e-006, is_small is true
my_test(0.001F) is true

my_test(0.001F) is false//] [/compare_floats_output]
*/
