
#include <iostream>
#include <iomanip>
#include "../performance/performance_measure.hpp"
#include "../performance/test_reference.cpp"

extern void polynomial_evaluate_0();
extern void polynomial_evaluate_1();
extern void polynomial_evaluate_2();
extern void polynomial_evaluate_3();
extern void polynomial_evaluate_m0();
extern void polynomial_evaluate_m1();
extern void polynomial_evaluate_m2();
extern void polynomial_evaluate_m3();
extern void rational_evaluate_0();
extern void rational_evaluate_1();
extern void rational_evaluate_2();
extern void rational_evaluate_3();
extern void rational_evaluate_m0();
extern void rational_evaluate_m1();
extern void rational_evaluate_m2();
extern void rational_evaluate_m3();

double total = 0;

int main(int argc)
{
   double reference_time = performance_measure(reference_evaluate);
   double test_time;

   std::cout << std::setprecision(3);

   //
   // BOOST_MATH_POLY_METHOD:
   //
   std::cout << std::setw(70) << std::left << std::scientific << "Testing BOOST_MATH_POLY_METHOD=0";
   test_time = performance_measure(polynomial_evaluate_0);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing BOOST_MATH_POLY_METHOD=1";
   test_time = performance_measure(polynomial_evaluate_1);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing BOOST_MATH_POLY_METHOD=2";
   test_time = performance_measure(polynomial_evaluate_2);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing BOOST_MATH_POLY_METHOD=3";
   test_time = performance_measure(polynomial_evaluate_3);
   std::cout << test_time - reference_time << std::endl;
   //
   // Mixed polynomial Arithmetic:
   //
   std::cout << std::setw(70) << std::left << std::scientific << "Testing int/double arithmetic with BOOST_MATH_POLY_METHOD=0";
   test_time = performance_measure(polynomial_evaluate_0);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing int/double arithmetic with BOOST_MATH_POLY_METHOD=1";
   test_time = performance_measure(polynomial_evaluate_1);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing int/double arithmetic with BOOST_MATH_POLY_METHOD=2";
   test_time = performance_measure(polynomial_evaluate_2);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing int/double arithmetic with BOOST_MATH_POLY_METHOD=3";
   test_time = performance_measure(polynomial_evaluate_3);
   std::cout << test_time - reference_time << std::endl;

   //
   // BOOST_MATH_RATIONAL_METHOD:
   //
   std::cout << std::setw(70) << std::left << std::scientific << "Testing BOOST_MATH_RATIONAL_METHOD=0";
   test_time = performance_measure(rational_evaluate_0);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing BOOST_MATH_RATIONAL_METHOD=1";
   test_time = performance_measure(rational_evaluate_1);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing BOOST_MATH_RATIONAL_METHOD=2";
   test_time = performance_measure(rational_evaluate_2);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing BOOST_MATH_RATIONAL_METHOD=3";
   test_time = performance_measure(rational_evaluate_3);
   std::cout << test_time - reference_time << std::endl;
   //
   // Mixed rational Arithmetic:
   //
   std::cout << std::setw(70) << std::left << std::scientific << "Testing int/double arithmetic with BOOST_MATH_RATIONAL_METHOD=0";
   test_time = performance_measure(rational_evaluate_0);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing int/double arithmetic with BOOST_MATH_RATIONAL_METHOD=1";
   test_time = performance_measure(rational_evaluate_1);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing int/double arithmetic with BOOST_MATH_RATIONAL_METHOD=2";
   test_time = performance_measure(rational_evaluate_2);
   std::cout << test_time - reference_time << std::endl;

   std::cout << std::setw(70) << std::left << std::scientific << "Testing int/double arithmetic with BOOST_MATH_RATIONAL_METHOD=3";
   test_time = performance_measure(rational_evaluate_3);
   std::cout << test_time - reference_time << std::endl;

   if(argc > 100000)
      std::cerr << total << std::endl;

   return 0;
}


void consume_result(double x)
{
   // Do nothing proceedure, don't let whole program optimisation
   // optimise this away - doing so may cause false readings....
   total += x;
}

void set_call_count(int i)
{}

