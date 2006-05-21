//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/ntl.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp> // for inverses
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/test.hpp>
#include <fstream>

#include <boost/math/tools/test_data.hpp>

#include "ntl_rr_lanczos.hpp"

using namespace boost::math::tools;
using namespace std;

struct erf_data_generator
{
   std::tr1::tuple<NTL::RR, NTL::RR> operator()(NTL::RR z)
   {
      // very naively calculate spots using the gamma function at high precision:
      int sign = 1;
      if(z < 0)
      {
         sign = -1;
         z = -z;
      }
      NTL::RR g1, g2;
      g1 = boost::math::tgamma_lower(NTL::RR(0.5), z * z);
      g1 /= sqrt(boost::math::constants::pi<NTL::RR>());
      g1 *= sign;

      if(z < 0.5)
      {
         g2 = 1 - (sign * g1);
      }
      else
      {
         g2 = boost::math::tgamma(NTL::RR(0.5), z * z);
         g2 /= sqrt(boost::math::constants::pi<NTL::RR>());
      }
      if(sign < 1)
         g2 = 2 - g2;
      return std::tr1::make_tuple(g1, g2);
   }
};

double double_factorial(int N)
{
   double result = 1;
   while(N > 2)
   {
      N -= 2;
      result *= N;
   }
   return result;
}

void asymptotic_limit(int Bits)
{
   //
   // The following block of code estimates how large z has
   // to be before we can use the asymptotic expansion for
   // erf/erfc and still get convergence: the series becomes
   // divergent eventually so we have to be careful!
   //
   double result = (std::numeric_limits<double>::max)();
   int terms = 0;
   for(int n = 1; n < 15; ++n)
   {
      double lim = (Bits-n) * log(2.0) - log(sqrt(3.14)) + log(double_factorial(2*n+1));
      double x = 1;
      while(x*x + (2*n+1)*log(x) <= lim)
         x += 0.1;
      if(x < result)
      {
         result = x;
         terms = n;
      }
   }

   std::cout << "Erf asymptotic limit for " 
      << Bits << " bit numbers is " 
      << result << " after approximately " 
      << terms << " terms." << std::endl;

   result = (std::numeric_limits<double>::max)();
   terms = 0;
   for(int n = 1; n < 30; ++n)
   {
      double x = pow(double_factorial(2*n+1)/pow(2.0, n-Bits), 1 / (2.0*n));
      if(x < result)
      {
         result = x;
         terms = n;
      }
   }

   std::cout << "Erfc asymptotic limit for " 
      << Bits << " bit numbers is " 
      << result << " after approximately " 
      << terms << " terms." << std::endl;
}

NTL::RR erfc_inv(NTL::RR r)
{
   std::cout << r << std::endl;
   return boost::math::erfc_inv(r);
}


int main(int argc, char* argv[])
{
   NTL::RR::SetPrecision(1000);
   NTL::RR::SetOutputPrecision(40);

   parameter_info<NTL::RR> arg1;
   test_data<NTL::RR> data;

   bool cont;
   std::string line;

   if(argc >= 2)
   {
      if(strcmp(argv[1], "--limits") == 0)
      {
         asymptotic_limit(24);
         asymptotic_limit(53);
         asymptotic_limit(64);
         asymptotic_limit(106);
         asymptotic_limit(113);
         return 0;
      }
      else if(strcmp(argv[1], "--erf_inv") == 0)
      {
         NTL::RR (*f)(NTL::RR);
         f = boost::math::erf_inv;
         std::cout << "Welcome.\n"
            "This program will generate spot tests for the inverse erf function:\n";
         std::cout << "Enter the number of data points: ";
         int points;
         std::cin >> points;
         data.insert(f, make_random_param(NTL::RR(-1), NTL::RR(1), points));
      }
      else if(strcmp(argv[1], "--erfc_inv") == 0)
      {
         NTL::RR (*f)(NTL::RR);
         f = /*boost::math::*/ erfc_inv;
         std::cout << "Welcome.\n"
            "This program will generate spot tests for the inverse erfc function:\n";
         std::cout << "Enter the number of data points: ";
         int points;
         std::cin >> points;
         data.insert(f, make_random_param(NTL::RR(0), NTL::RR(2), points));
      }
   }
   else
   {
      std::cout << "Welcome.\n"
         "This program will generate spot tests for the erf and erfc functions:\n"
         "  erf(z) and erfc(z)\n\n";

      do{
         if(0 == get_user_parameter_info(arg1, "a"))
            return 1;
         data.insert(erf_data_generator(), arg1);

         std::cout << "Any more data [y/n]?";
         std::getline(std::cin, line);
         boost::algorithm::trim(line);
         cont = (line == "y");
      }while(cont);
   }

   std::cout << "Enter name of test data file [default=erf_data.ipp]";
   std::getline(std::cin, line);
   boost::algorithm::trim(line);
   if(line == "")
      line = "erf_data.ipp";
   std::ofstream ofs(line.c_str());
   write_code(ofs, data, "erf_data");
   
   return 0;
}

/* Output for asymptotic limits:

Erf asymptotic limit for 24 bit numbers is 2.8 after approximately 6 terms.
Erfc asymptotic limit for 24 bit numbers is 4.12064 after approximately 17 terms.
Erf asymptotic limit for 53 bit numbers is 4.3 after approximately 11 terms.
Erfc asymptotic limit for 53 bit numbers is 6.19035 after approximately 29 terms.
Erf asymptotic limit for 64 bit numbers is 4.8 after approximately 12 terms.
Erfc asymptotic limit for 64 bit numbers is 7.06004 after approximately 29 terms.
Erf asymptotic limit for 106 bit numbers is 6.5 after approximately 14 terms.
Erfc asymptotic limit for 106 bit numbers is 11.6626 after approximately 29 terms.
Erf asymptotic limit for 113 bit numbers is 6.8 after approximately 14 terms.
Erfc asymptotic limit for 113 bit numbers is 12.6802 after approximately 29 terms.
*/

