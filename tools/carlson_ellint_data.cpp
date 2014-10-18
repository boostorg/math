// Copyright John Maddock 2006.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/test_data.hpp>
#include <boost/test/included/prg_exec_monitor.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <fstream>
#include <boost/math/tools/test_data.hpp>
#include <boost/random.hpp>
#include "mp_t.hpp"

float extern_val;
// confuse the compilers optimiser, and force a truncation to float precision:
float truncate_to_float(float const * pf)
{
   extern_val = *pf;
   return *pf;
}

boost::math::tuple<mp_t, mp_t, mp_t, mp_t> generate_rf_data(mp_t n)
{
   static boost::mt19937 r;
   boost::uniform_real<float> ur(0, 1);
   boost::uniform_int<int> ui(-100, 100);
   float x = ur(r);
   x = ldexp(x, ui(r));
   mp_t xr(truncate_to_float(&x));
   float y = ur(r);
   y = ldexp(y, ui(r));
   mp_t yr(truncate_to_float(&y));
   float z = ur(r);
   z = ldexp(z, ui(r));
   mp_t zr(truncate_to_float(&z));

   mp_t result = boost::math::ellint_rf(xr, yr, zr);
   return boost::math::make_tuple(xr, yr, zr, result);
}

boost::math::tuple<mp_t, mp_t, mp_t> generate_rc_data(mp_t n)
{
   static boost::mt19937 r;
   boost::uniform_real<float> ur(0, 1);
   boost::uniform_int<int> ui(-100, 100);
   float x = ur(r);
   x = ldexp(x, ui(r));
   mp_t xr(truncate_to_float(&x));
   float y = ur(r);
   y = ldexp(y, ui(r));
   mp_t yr(truncate_to_float(&y));

   mp_t result = boost::math::ellint_rc(xr, yr);
   return boost::math::make_tuple(xr, yr, result);
}

boost::math::tuple<mp_t, mp_t, mp_t, mp_t, mp_t> generate_rj_data(mp_t n)
{
   static boost::mt19937 r;
   boost::uniform_real<float> ur(0, 1);
   boost::uniform_real<float> nur(-1, 1);
   boost::uniform_int<int> ui(-100, 100);
   float x = ur(r);
   x = ldexp(x, ui(r));
   mp_t xr(truncate_to_float(&x));
   float y = ur(r);
   y = ldexp(y, ui(r));
   mp_t yr(truncate_to_float(&y));
   float z = ur(r);
   z = ldexp(z, ui(r));
   mp_t zr(truncate_to_float(&z));
   float p = nur(r);
   p = ldexp(p, ui(r));
   mp_t pr(truncate_to_float(&p));

   boost::math::ellint_rj(x, y, z, p);

   mp_t result = boost::math::ellint_rj(xr, yr, zr, pr);
   return boost::math::make_tuple(xr, yr, zr, pr, result);
}

boost::math::tuple<mp_t, mp_t, mp_t, mp_t> generate_rd_data(mp_t n)
{
   static boost::mt19937 r;
   boost::uniform_real<float> ur(0, 1);
   boost::uniform_int<int> ui(-100, 100);
   float x = ur(r);
   x = ldexp(x, ui(r));
   mp_t xr(truncate_to_float(&x));
   float y = ur(r);
   y = ldexp(y, ui(r));
   mp_t yr(truncate_to_float(&y));
   float z = ur(r);
   z = ldexp(z, ui(r));
   mp_t zr(truncate_to_float(&z));

   mp_t result = boost::math::ellint_rd(xr, yr, zr);
   return boost::math::make_tuple(xr, yr, zr, result);
}

int cpp_main(int argc, char*argv [])
{
   using namespace boost::math::tools;

   parameter_info<mp_t> arg1, arg2;
   test_data<mp_t> data;

   bool cont;
   std::string line;

   if(argc < 1)
      return 1;

   do{
      int count;
      std::cout << "Number of points: ";
      std::cin >> count;
      
      arg1 = make_periodic_param(mp_t(0), mp_t(1), count);
      arg1.type |= dummy_param;

      //
      // Change this next line to get the R variant you want:
      //
      data.insert(&generate_rd_data, arg1);

      std::cout << "Any more data [y/n]?";
      std::getline(std::cin, line);
      boost::algorithm::trim(line);
      cont = (line == "y");
   }while(cont);

   std::cout << "Enter name of test data file [default=ellint_rf_data.ipp]";
   std::getline(std::cin, line);
   boost::algorithm::trim(line);
   if(line == "")
      line = "ellint_rf_data.ipp";
   std::ofstream ofs(line.c_str());
   line.erase(line.find('.'));
   ofs << std::scientific << std::setprecision(40);
   write_code(ofs, data, line.c_str());

   return 0;
}


