
#include <boost/math/tools/ntl.hpp>
#include <boost/limits.hpp>
#include <vector>

void write_table(unsigned max_exponent)
{
   NTL::RR max = ldexp(NTL::RR(1), max_exponent);

   std::vector<NTL::RR> factorials;
   factorials.push_back(1);

   NTL::RR f(1);
   unsigned i = 1;

   while(f < max)
   {
      factorials.push_back(f);
      ++i;
      f *= i;
   }

   //
   // now write out the results to cout:
   //
   std::cout << std::scientific;
   std::cout << "   static const boost::array<T, " << factorials.size() << "> factorials = {\n";
   for(unsigned j = 0; j < factorials.size(); ++j)
      std::cout << "      " << factorials[j] << "L,\n";
   std::cout << "   };\n\n";
}


int main()
{
   NTL::RR::SetPrecision(300);
   NTL::RR::SetOutputPrecision(40);
   write_table(16384/*std::numeric_limits<float>::max_exponent*/);
}
