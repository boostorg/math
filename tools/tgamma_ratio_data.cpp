//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/ntl.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/test.hpp>
#include <fstream>

#include <boost/math/tools/test_data.hpp>

#include "ntl_rr_lanczos.hpp"

using namespace boost::math::tools;
using namespace std;

std::tr1::tuple<NTL::RR, NTL::RR> 
   tgamma_ratio(const NTL::RR& a, const NTL::RR& delta)
{
   if(delta > a)
      throw std::domain_error("");
   NTL::RR tg = boost::math::tgamma(a);
   NTL::RR r1 = tg / boost::math::tgamma(a + delta);
   NTL::RR r2 = tg / boost::math::tgamma(a - delta);
   if((r1 > (std::numeric_limits<float>::max)()) || (r2 > (std::numeric_limits<float>::max)()))
      throw std::domain_error("");

   return std::tr1::make_tuple(r1, r2);
}

NTL::RR tgamma_ratio2(const NTL::RR& a, const NTL::RR& b)
{
   return boost::math::tgamma(a) / boost::math::tgamma(b);
}


int test_main(int argc, char*argv [])
{
   NTL::RR::SetPrecision(1000);
   NTL::RR::SetOutputPrecision(40);

   parameter_info<NTL::RR> arg1, arg2;
   test_data<NTL::RR> data;

   bool cont;
   std::string line;

   if((argc >= 2) && (strcmp(argv[1], "--ratio") == 0))
   {
      std::cout << "Welcome.\n"
         "This program will generate spot tests for the function tgamma_ratio(a, b)\n\n";

      do{
         if(0 == get_user_parameter_info(arg1, "a"))
            return 1;
         if(0 == get_user_parameter_info(arg2, "b"))
            return 1;
         data.insert(&tgamma_ratio2, arg1, arg2);

         std::cout << "Any more data [y/n]?";
         std::getline(std::cin, line);
         boost::algorithm::trim(line);
         cont = (line == "y");
      }while(cont);
   }
   else
   {
      std::cout << "Welcome.\n"
         "This program will generate spot tests for the function tgamma_delta_ratio(a, delta)\n\n";

      do{
         if(0 == get_user_parameter_info(arg1, "a"))
            return 1;
         if(0 == get_user_parameter_info(arg2, "delta"))
            return 1;
         data.insert(&tgamma_ratio, arg1, arg2);

         std::cout << "Any more data [y/n]?";
         std::getline(std::cin, line);
         boost::algorithm::trim(line);
         cont = (line == "y");
      }while(cont);
   }

   std::cout << "Enter name of test data file [default=tgamma_ratio_data.ipp]";
   std::getline(std::cin, line);
   boost::algorithm::trim(line);
   if(line == "")
      line = "tgamma_ratio_data.ipp";
   std::ofstream ofs(line.c_str());
   ofs << std::scientific;
   write_code(ofs, data, "tgamma_ratio_data");
   
   return 0;
}


