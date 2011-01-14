//  Copyright John Maddock 2007.
//  Copyright Paul A. Bristow 2010
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Note that this file contains quickbook mark-up as well as code
// and comments, don't change any of the special comment mark-ups!

//[policy_eg_10

/*`

To understand how the rounding policies for 
the discrete distributions can be used, we'll
use the 50-sample binomial distribution with a 
success fraction of 0.5 once again, and calculate
all the possible quantiles at 0.05 and 0.95.

Begin by including the needed headers (and some using statements for conciseness):

*/
#include <ios>
#include <iostream>
using std::cout; using std::endl;
using std::left; using std::fixed; using std::right; using std::scientific;
#include <boost/detail/iomanip.hpp>
using boost::detail::setw;
using boost::detail::setprecision;

#include <boost/math/distributions/binomial.hpp>
/*`

Next we'll bring the needed declarations into scope, and
define distribution types for all the available rounding policies:

*/
// Avoid 
// using namespace std; // and 
// using namespace boost::math;
// to avoid potential ambiguity of names, like binomial.
// using namespace boost::math::policies; is small risk, but
// the necessary items are brought into scope thus:

using boost::math::binomial_distribution;
using boost::math::policies::policy;
using boost::math::policies::discrete_quantile;

using boost::math::policies::integer_round_outwards;
using boost::math::policies::integer_round_down;
using boost::math::policies::integer_round_up;
using boost::math::policies::integer_round_nearest;
using boost::math::policies::integer_round_inwards;
using boost::math::policies::real;

using boost::math::binomial_distribution; // Not std::binomial_distribution.

typedef binomial_distribution<
            double, 
            policy<discrete_quantile<integer_round_outwards> > > 
        binom_round_outwards;

typedef binomial_distribution<
            double, 
            policy<discrete_quantile<integer_round_inwards> > > 
        binom_round_inwards;

typedef binomial_distribution<
            double, 
            policy<discrete_quantile<integer_round_down> > > 
        binom_round_down;

typedef binomial_distribution<
            double, 
            policy<discrete_quantile<integer_round_up> > > 
        binom_round_up;

typedef binomial_distribution<
            double, 
            policy<discrete_quantile<integer_round_nearest> > > 
        binom_round_nearest;

typedef binomial_distribution<
            double, 
            policy<discrete_quantile<real> > > 
        binom_real_quantile;

/*`
Now let's set to work calling those quantiles:
*/

int main()
{
   cout << 
      "Testing rounding policies for a 50 sample binomial distribution,\n"
      "with a success fraction of 0.5.\n\n"
      "Lower quantiles are calculated at p = 0.05\n\n"
      "Upper quantiles at p = 0.95.\n\n";

   cout << boost::detail::setw(25) << right
      << "Policy"<< boost::detail::setw(18) << right 
      << "Lower Quantile" << boost::detail::setw(18) << right 
      << "Upper Quantile" << endl;
   
   // Test integer_round_outwards:
   cout << boost::detail::setw(25) << right
      << "integer_round_outwards"
      << boost::detail::setw(18) << right
      << quantile(binom_round_outwards(50, 0.5), 0.05)
      << boost::detail::setw(18) << right
      << quantile(binom_round_outwards(50, 0.5), 0.95) 
      << endl;
   
   // Test integer_round_inwards:
   cout << boost::detail::setw(25) << right
      << "integer_round_inwards"
      << boost::detail::setw(18) << right
      << quantile(binom_round_inwards(50, 0.5), 0.05)
      << boost::detail::setw(18) << right
      << quantile(binom_round_inwards(50, 0.5), 0.95) 
      << endl;
   
   // Test integer_round_down:
   cout << boost::detail::setw(25) << right
      << "integer_round_down"
      << boost::detail::setw(18) << right
      << quantile(binom_round_down(50, 0.5), 0.05)
      << boost::detail::setw(18) << right
      << quantile(binom_round_down(50, 0.5), 0.95) 
      << endl;
   
   // Test integer_round_up:
   cout << boost::detail::setw(25) << right
      << "integer_round_up"
      << boost::detail::setw(18) << right
      << quantile(binom_round_up(50, 0.5), 0.05)
      << boost::detail::setw(18) << right
      << quantile(binom_round_up(50, 0.5), 0.95) 
      << endl;
   
   // Test integer_round_nearest:
   cout << boost::detail::setw(25) << right
      << "integer_round_nearest"
      << boost::detail::setw(18) << right
      << quantile(binom_round_nearest(50, 0.5), 0.05)
      << boost::detail::setw(18) << right
      << quantile(binom_round_nearest(50, 0.5), 0.95) 
      << endl;
   
   // Test real:
   cout << boost::detail::setw(25) << right
      << "real"
      << boost::detail::setw(18) << right
      << quantile(binom_real_quantile(50, 0.5), 0.05)
      << boost::detail::setw(18) << right
      << quantile(binom_real_quantile(50, 0.5), 0.95) 
      << endl;
} // int main()

/*`

Which produces the program output:

[pre
  policy_eg_10.vcxproj -> J:\Cpp\MathToolkit\test\Math_test\Release\policy_eg_10.exe
  Testing rounding policies for a 50 sample binomial distribution,
  with a success fraction of 0.5.
  
  Lower quantiles are calculated at p = 0.05
  
  Upper quantiles at p = 0.95.
  
                     Policy    Lower Quantile    Upper Quantile
     integer_round_outwards                18                31
      integer_round_inwards                19                30
         integer_round_down                18                30
           integer_round_up                19                31
      integer_round_nearest                19                30
                       real            18.701            30.299
]

*/

//] //[policy_eg_10] ends quickbook import.
