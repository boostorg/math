
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include "performance_measure.hpp"

extern void reference_evaluate();
/*
extern void polynomial_evaluate();
extern void polynomial_mixed_evaluate();
extern void rational_evaluate();
extern void rational_mixed_evaluate();
extern void gamma_evaluate();
extern void lgamma_evaluate();
extern void erf_evaluate();
extern void igamma_evaluate();
extern void igamma_inv_evaluate();
extern void ibeta_evaluate();
extern void ibeta_inv_evaluate();



test_info info[] = {
   { "polynomial", &polynomial_evaluate },
   { "polynomial_mixed", &polynomial_mixed_evaluate },
   { "rational", &rational_evaluate },
   { "rational_mixed", &rational_mixed_evaluate },
   { "gamma", &gamma_evaluate },
   { "lgamma", &lgamma_evaluate },
   { "erf", &erf_evaluate },
   { "igamma_inv", &igamma_inv_evaluate },
   { "igamma", &igamma_evaluate },
   { "ibeta_inv", &ibeta_inv_evaluate },
   { "ibeta", &ibeta_evaluate },
};
*/

std::map<std::string, double> times;
std::set<test_info> tests;
double total = 0;
int call_count = 0;

std::set<test_info>& all_tests()
{
   static std::set<test_info> i;
   return i;
}

void add_new_test(test_info i)
{
   all_tests().insert(i);
}

void set_call_count(int i)
{
   call_count = i;
}

void show_help()
{
   std::cout << "Specify on the command line which functions to test.\n"
      "Available options are:\n";
   std::set<test_info>::const_iterator i(all_tests().begin()), j(all_tests().end());
   while(i != j)
   {
      std::cout << "   --" << (*i).name << std::endl;
      ++i;
   }
   std::cout << "Or use --all to test everything." << std::endl;
}

void run_tests()
{
   // Get time for empty proceedure:
   double reference_time = performance_measure(reference_evaluate);

   std::set<test_info>::const_iterator i, j;
   for(i = tests.begin(), j = tests.end(); i != j; ++i)
   {
      set_call_count(1);
      std::cout << "Testing " << std::left << std::setw(40) << i->name << std::flush;
      double time = performance_measure(i->proc) - reference_time;
      time /= call_count;
      std::cout << std::setprecision(3) << std::scientific << time << std::endl;
   }
}

int main(int argc, const char** argv)
{
   try{

   if(argc >= 2)
   {
      for(int i = 1; i < argc; ++i)
      {
         if(std::strcmp(argv[i], "--help") == 0)
         {
            show_help();
         }
         else if(std::strcmp(argv[i], "--all") == 0)
         {
            std::set<test_info>::const_iterator a(all_tests().begin()), b(all_tests().end());
            while(a != b)
            {
               tests.insert(*a);
               ++a;
            }
         }
         else
         {
            bool found = false;
            if((argv[i][0] == '-') && (argv[i][1] == '-'))
            {
               std::set<test_info>::const_iterator a(all_tests().begin()), b(all_tests().end());
               while(a != b)
               {
                  if(std::strcmp(argv[i] + 2, (*a).name) == 0)
                  {
                     found = true;
                     tests.insert(*a);
                     break;
                  }
                  ++a;
               }
            }
            if(!found)
            {
               std::cerr << "Unknown option: " << argv[i] << std::endl;
               return 1;
            }
         }
      }
   }
   else
   {
      show_help();
   }
   run_tests();
   //
   // This is just to confuse the optimisers:
   //
   if(argc == 100000)
   {
      std::cerr << total << std::endl;
   }

   }
   catch(const std::exception& e)
   {
      std::cerr << e.what() << std::endl;
   }

   return 0;
}

void consume_result(double x)
{
   // Do nothing proceedure, don't let whole program optimisation
   // optimise this away - doing so may cause false readings....
   total += x;
}
