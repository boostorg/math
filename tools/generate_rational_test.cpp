
#include <boost/math/tools/ntl.hpp>
#include <boost/tr1/random.hpp>
#include <boost/math/tools/rational.hpp>
#include <iostream>
#include <fstream>

int main()
{
   using namespace boost::math;
   using namespace boost::math::tools;

   NTL::RR::SetPrecision(500);
   NTL::RR::SetOutputPrecision(40);

   std::tr1::mt19937 rnd;
   std::tr1::variate_generator<
      std::tr1::mt19937, 
      std::tr1::uniform_int<> > gen(rnd, std::tr1::uniform_int<>(1, 12));

   for(unsigned i = 1; i < 12; ++i)
   {
      std::vector<int> coef;
      for(unsigned j = 0; j < i; ++j)
      {
         coef.push_back(gen());
      }
      std::cout << 
"   //\n"
"   // Polynomials of order " << i-1 << "\n"
"   //\n"
"   static const U n" << i << "c[" << i << "] = { ";
      for(unsigned j = 0; j < i; ++j)
      {
         if(j) 
            std::cout  << ", ";
         std::cout << coef[j];
      }
      std::cout << " };\n";
      std::cout <<
         "   static const boost::array<U, " << i << "> n" << i << "a = { ";
      for(unsigned j = 0; j < i; ++j)
      {
         if(j) 
            std::cout  << ", ";
         std::cout << coef[j];
      }
      std::cout << " };\n";

      NTL::RR r1 = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(0.125), i);
      NTL::RR r2 = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(0.25), i);
      NTL::RR r3 = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(0.75), i);
      NTL::RR r4 = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(1) - NTL::RR(1) / 64, i);

      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "c, static_cast<T>(0.125), " << i << "),\n"
         "      static_cast<T>(" << r1 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "c, static_cast<T>(0.25), " << i << "),\n"
         "      static_cast<T>(" << r2 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "c, static_cast<T>(0.75), " << i << "),\n"
         "      static_cast<T>(" << r3 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "c, static_cast<T>(1.0f - 1.0f/64.0f), " << i << "),\n"
         "      static_cast<T>(" << r4 << "L),\n"
         "      tolerance);\n\n";

      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "c, static_cast<T>(0.125)),\n"
         "      static_cast<T>(" << r1 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "c, static_cast<T>(0.25)),\n"
         "      static_cast<T>(" << r2 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "c, static_cast<T>(0.75)),\n"
         "      static_cast<T>(" << r3 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "c, static_cast<T>(1.0f - 1.0f/64.0f)),\n"
         "      static_cast<T>(" << r4 << "L),\n"
         "      tolerance);\n\n";

      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "a, static_cast<T>(0.125)),\n"
         "      static_cast<T>(" << r1 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "a, static_cast<T>(0.25)),\n"
         "      static_cast<T>(" << r2 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "a, static_cast<T>(0.75)),\n"
         "      static_cast<T>(" << r3 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_polynomial(n" << i << "a, static_cast<T>(1.0f - 1.0f/64.0f)),\n"
         "      static_cast<T>(" << r4 << "L),\n"
         "      tolerance);\n\n";

      r1 = boost::math::tools::evaluate_even_polynomial(&coef[0], NTL::RR(0.125), i);
      r2 = boost::math::tools::evaluate_even_polynomial(&coef[0], NTL::RR(0.25), i);
      r3 = boost::math::tools::evaluate_even_polynomial(&coef[0], NTL::RR(0.75), i);
      r4 = boost::math::tools::evaluate_even_polynomial(&coef[0], NTL::RR(1) - NTL::RR(1) / 64, i);

      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "c, static_cast<T>(0.125), " << i << "),\n"
         "      static_cast<T>(" << r1 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "c, static_cast<T>(0.25), " << i << "),\n"
         "      static_cast<T>(" << r2 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "c, static_cast<T>(0.75), " << i << "),\n"
         "      static_cast<T>(" << r3 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "c, static_cast<T>(1.0f - 1.0f/64.0f), " << i << "),\n"
         "      static_cast<T>(" << r4 << "L),\n"
         "      tolerance);\n\n";

      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "c, static_cast<T>(0.125)),\n"
         "      static_cast<T>(" << r1 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "c, static_cast<T>(0.25)),\n"
         "      static_cast<T>(" << r2 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "c, static_cast<T>(0.75)),\n"
         "      static_cast<T>(" << r3 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "c, static_cast<T>(1.0f - 1.0f/64.0f)),\n"
         "      static_cast<T>(" << r4 << "L),\n"
         "      tolerance);\n\n";

      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "a, static_cast<T>(0.125)),\n"
         "      static_cast<T>(" << r1 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "a, static_cast<T>(0.25)),\n"
         "      static_cast<T>(" << r2 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "a, static_cast<T>(0.75)),\n"
         "      static_cast<T>(" << r3 << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_even_polynomial(n" << i << "a, static_cast<T>(1.0f - 1.0f/64.0f)),\n"
         "      static_cast<T>(" << r4 << "L),\n"
         "      tolerance);\n\n";

      if(i > 1)
      {
         r1 = boost::math::tools::evaluate_odd_polynomial(&coef[0], NTL::RR(0.125), i);
         r2 = boost::math::tools::evaluate_odd_polynomial(&coef[0], NTL::RR(0.25), i);
         r3 = boost::math::tools::evaluate_odd_polynomial(&coef[0], NTL::RR(0.75), i);
         r4 = boost::math::tools::evaluate_odd_polynomial(&coef[0], NTL::RR(1) - NTL::RR(1) / 64, i);

         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "c, static_cast<T>(0.125), " << i << "),\n"
            "      static_cast<T>(" << r1 << "L),\n"
            "      tolerance);\n";
         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "c, static_cast<T>(0.25), " << i << "),\n"
            "      static_cast<T>(" << r2 << "L),\n"
            "      tolerance);\n";
         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "c, static_cast<T>(0.75), " << i << "),\n"
            "      static_cast<T>(" << r3 << "L),\n"
            "      tolerance);\n";
         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "c, static_cast<T>(1.0f - 1.0f/64.0f), " << i << "),\n"
            "      static_cast<T>(" << r4 << "L),\n"
            "      tolerance);\n\n";

         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "c, static_cast<T>(0.125)),\n"
            "      static_cast<T>(" << r1 << "L),\n"
            "      tolerance);\n";
         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "c, static_cast<T>(0.25)),\n"
            "      static_cast<T>(" << r2 << "L),\n"
            "      tolerance);\n";
         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "c, static_cast<T>(0.75)),\n"
            "      static_cast<T>(" << r3 << "L),\n"
            "      tolerance);\n";
         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "c, static_cast<T>(1.0f - 1.0f/64.0f)),\n"
            "      static_cast<T>(" << r4 << "L),\n"
            "      tolerance);\n\n";

         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "a, static_cast<T>(0.125)),\n"
            "      static_cast<T>(" << r1 << "L),\n"
            "      tolerance);\n";
         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "a, static_cast<T>(0.25)),\n"
            "      static_cast<T>(" << r2 << "L),\n"
            "      tolerance);\n";
         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "a, static_cast<T>(0.75)),\n"
            "      static_cast<T>(" << r3 << "L),\n"
            "      tolerance);\n";
         std::cout <<
            "   BOOST_CHECK_CLOSE(\n"
            "      boost::math::tools::evaluate_odd_polynomial(n" << i << "a, static_cast<T>(1.0f - 1.0f/64.0f)),\n"
            "      static_cast<T>(" << r4 << "L),\n"
            "      tolerance);\n\n";
      }

      r1 = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(0.125), i);
      r2 = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(0.25), i);
      r3 = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(0.75), i);
      r4 = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(1) - NTL::RR(1) / 64, i);

      coef.clear();
      for(unsigned j = 0; j < i; ++j)
      {
         coef.push_back(gen());
      }
      std::cout << 
"   //\n"
"   // Rational functions of order " << i-1 << "\n"
"   //\n"
"   static const U d" << i << "c[" << i << "] = { ";
      for(unsigned j = 0; j < i; ++j)
      {
         if(j) 
            std::cout  << ", ";
         std::cout << coef[j];
      }
      std::cout << " };\n";
      std::cout <<
         "   static const boost::array<U, " << i << "> d" << i << "a = { ";
      for(unsigned j = 0; j < i; ++j)
      {
         if(j) 
            std::cout  << ", ";
         std::cout << coef[j];
      }
      std::cout << " };\n";

      NTL::RR r1d = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(0.125), i);
      NTL::RR r2d = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(0.25), i);
      NTL::RR r3d = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(0.75), i);
      NTL::RR r4d = boost::math::tools::evaluate_polynomial(&coef[0], NTL::RR(1) - NTL::RR(1) / 64, i);

      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "c, d" << i << "c, static_cast<T>(0.125), " << i << "),\n"
         "      static_cast<T>(" << r1/r1d << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "c, d" << i << "c, static_cast<T>(0.25), " << i << "),\n"
         "      static_cast<T>(" << r2/r2d << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "c, d" << i << "c, static_cast<T>(0.75), " << i << "),\n"
         "      static_cast<T>(" << r3/r3d << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "c, d" << i << "c, static_cast<T>(1.0f - 1.0f/64.0f), " << i << "),\n"
         "      static_cast<T>(" << r4/r4d << "L),\n"
         "      tolerance);\n\n";

      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "c, d" << i << "c, static_cast<T>(0.125)),\n"
         "      static_cast<T>(" << r1/r1d << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "c, d" << i << "c, static_cast<T>(0.25)),\n"
         "      static_cast<T>(" << r2/r2d << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "c, d" << i << "c, static_cast<T>(0.75)),\n"
         "      static_cast<T>(" << r3/r3d << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "c, d" << i << "c, static_cast<T>(1.0f - 1.0f/64.0f)),\n"
         "      static_cast<T>(" << r4/r4d << "L),\n"
         "      tolerance);\n\n";

      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "a, d" << i << "a, static_cast<T>(0.125)),\n"
         "      static_cast<T>(" << r1/r1d << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "a, d" << i << "a, static_cast<T>(0.25)),\n"
         "      static_cast<T>(" << r2/r2d << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "a, d" << i << "a, static_cast<T>(0.75)),\n"
         "      static_cast<T>(" << r3/r3d << "L),\n"
         "      tolerance);\n";
      std::cout <<
         "   BOOST_CHECK_CLOSE(\n"
         "      boost::math::tools::evaluate_rational(n" << i << "a, d" << i << "a, static_cast<T>(1.0f - 1.0f/64.0f)),\n"
         "      static_cast<T>(" << r4/r4d << "L),\n"
         "      tolerance);\n\n";
   }

   return 0;
}



