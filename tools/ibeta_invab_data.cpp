//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/ntl.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>

#include <boost/math/tools/test_data.hpp>
#include "ntl_rr_lanczos.hpp"

using namespace boost::math::tools;

//
// Force trunctation to float precision of input values:
// we must ensure that the input values are exactly representable
// in whatever type we are testing, or the output values will all
// be thrown off:
//
float external_f;
float force_truncate(const float* f)
{
   external_f = *f;
   return external_f;
}

float truncate_to_float(NTL::RR r)
{
   float f = boost::math::tools::real_cast<float>(r);
   return force_truncate(&f);
}

std::tr1::mt19937 rnd;
std::tr1::uniform_real<float> ur_a(1.0F, 5.0F);
std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<float> > gen(rnd, ur_a);
std::tr1::uniform_real<float> ur_a2(0.0F, 100.0F);
std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<float> > gen2(rnd, ur_a2);

struct ibeta_inv_data_generator
{
   std::tr1::tuple<NTL::RR, NTL::RR, NTL::RR, NTL::RR, NTL::RR, NTL::RR, NTL::RR> operator()
      (NTL::RR bp, NTL::RR x_, NTL::RR p_)
   {
      float b = truncate_to_float(real_cast<float>(gen() * pow(NTL::RR(10), bp))); 
      float x = truncate_to_float(real_cast<float>(x_));
      float p = truncate_to_float(real_cast<float>(p_));
      std::cout << b << " " << x << " " << p << std::flush;
      NTL::RR inv = boost::math::ibeta_inva(NTL::RR(b), NTL::RR(x), NTL::RR(p));
      std::cout << " " << inv << std::flush;
      NTL::RR invc = boost::math::ibetac_inva(NTL::RR(b), NTL::RR(x), NTL::RR(p));
      std::cout << " " << invc << std::endl;
      NTL::RR invb = boost::math::ibeta_invb(NTL::RR(b), NTL::RR(x), NTL::RR(p));
      std::cout << " " << invb << std::flush;
      NTL::RR invbc = boost::math::ibetac_invb(NTL::RR(b), NTL::RR(x), NTL::RR(p));
      std::cout << " " << invbc << std::endl;
      return std::tr1::make_tuple(b, x, p, inv, invc, invb, invbc);
   }
};

int test_main(int argc, char*argv [])
{
   NTL::RR::SetPrecision(1000);
   NTL::RR::SetOutputPrecision(100);

   bool cont;
   std::string line;

   parameter_info<NTL::RR> arg1, arg2, arg3;
   test_data<NTL::RR> data;

   std::cout << "Welcome.\n"
      "This program will generate spot tests for the inverse incomplete beta function:\n"
      "  ibeta_inva(a, p) and ibetac_inva(a, q)\n\n";

   arg1 = make_periodic_param(NTL::RR(-5), NTL::RR(6), 11);
   arg2 = make_random_param(NTL::RR(0.0001), NTL::RR(1), 10);
   arg3 = make_random_param(NTL::RR(0.0001), NTL::RR(1), 10);

   arg1.type |= dummy_param;
   arg2.type |= dummy_param;
   arg3.type |= dummy_param;

   data.insert(ibeta_inv_data_generator(), arg1, arg2, arg3);

   line = "ibeta_inva_data.ipp";
   std::ofstream ofs(line.c_str());
   write_code(ofs, data, "ibeta_inva_data");
   
   return 0;
}

