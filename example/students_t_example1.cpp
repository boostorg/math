// students_t_example1.cpp
// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Example 1 of using Student's t

// http://en.wikipedia.org/wiki/Student's_t-test  says:
// The t statistic was invented by William Sealy Gosset
// for cheaply monitoring the quality of beer brews.
// "Student" was his pen name.
// WS Gosset was statistician for Guinness brewery in Dublin, Ireland,
// hired due to Claude Guinness's innovative policy of recruiting the
// best graduates from Oxford and Cambridge for applying biochemistry
// and statistics to Guinness's industrial processes.
// Gosset published the t test in Biometrika in 1908,
// but was forced to use a pen name by his employer who regarded the fact
// that they were using statistics as a trade secret.
// In fact, Gosset's identity was unknown not only to fellow statisticians
// but to his employer - the company insisted on the pseudonym
// so that it could turn a blind eye to the breach of its rules.

// Data for this example from:
// P.K.Hou, O. W. Lau & M.C. Wong, Analyst (1983) vol. 108, p 64.
// from Statistics for Analytical Chemistry, 3rd ed. (1994), pp 54-55
// J. C. Miller and J. N. Miller, Ellis Horwood ISBN 0 13 0309907

// Determination of mercury by cold-vapour atomic absorption,
// the following values were obtained fusing a trusted
// Standard Reference Material containing 38.9% mercury,
// which we assume is correct or 'true'.
double standard = 38.9;

const int values = 3;
double value[values] = {38.9, 37.4, 37.1};

// Is there any evidence for systematic error?

// http://en.wikipedia.org/wiki/Student%27s_t_distribution

#include <boost/math/special_functions/students_t.hpp>
	using boost::math::students_t;  // Probability of students_t(df, t).

#include <iostream>
	using std::cout;
	using std::endl;
#include <iomanip>
	using std::setprecision;
#include <cmath>
	using std::sqrt;

int main()
{
	cout << "Example 1 using Student's t function. ";
#if defined(__FILE__) && defined(__TIMESTAMP__)
	cout << "  " << __FILE__ << ' ' << __TIMESTAMP__ << ' '<< _MSC_FULL_VER;
#endif
	cout << endl;

	// Example/test using tabulated value (deliberately coded as naively as possible).

	// Null hypothesis is that there is no difference (greater or less) between measured and standard.

	double degrees_of_freedom = values-1; // 3-1 = 2

	double mean = (value[0] + value[1] + value[2]) / static_cast<double>(values);
	cout << "mean = " << mean << ", (mean - standard) = " << mean - standard  << endl;
	double sd = sqrt(((value[0] - mean) * (value[0] - mean) + (value[1] - mean) * (value[1] - mean) + (value[2] - mean) * (value[2] - mean))/ static_cast<double>(values-1));
	cout << "standard deviation = " << sd << endl;
	if (sd == 0.)
	{
		cout << "Measured mean is identical to SRM value," << endl;
		cout << "so probability of no difference between measured and standard (the 'null hypothesis') is unity." << endl;
		return 0;
	}
	
	double t = (mean - standard) * std::sqrt(static_cast<double>(values)) / sd;
	cout << "Student's t = " << t << endl;

	cout.precision(5); // Useful accuracy is only a few decimal digits, but seems to give at least 5.
	cout << "Probability of Student's t is " << students_t(degrees_of_freedom, abs(t)) << endl; //  0.90657, is 1 tailed.
	// So there is insufficient evidence of a difference to meet a 95% (1 in 20) criterion.

	return 0;
}  // int main()

/*

Output is:
Example 1 using Student's t function.   ..\..\..\..\..\..\boost-sandbox\libs\math_functions\example\students_t_example1.cpp Sun Jun 18 12:37:42 2006 140050727
mean = 37.8, (mean - standard) = -1.1
standard deviation = 0.964365
Student's t = -1.97566
Probability of Student's t is 0.90657
Press any key to continue . . .


*/

/*
Obselete junk?
	double p = students_t(degrees_of_freedom, t);
	p = students_t(degrees_of_freedom, t);
	cout << "Probability of t = " << t << " is " << p << endl;

	//t = 4.30; 4.30265272974946
	//p = students_t(degrees_of_freedom, t); // 0.975 = tquantiles value
	//cout << "Probability of t = " << t << " is " << p << endl;
	// But MathCAD pt(2, 4.30265272974946) = 0.944426429950027 ??? expect 0.95

	//t = 2.92;
	//p = students_t(degrees_of_freedom, t); // 0.95
	//cout << "Probability of t = " << t << " is " << p << endl;
	//t = -2.92;
	//p = students_t(degrees_of_freedom, t); // 0.05
	//cout << "Probability of t = " << t << " is " << p << endl;

	// Means probability(t < 1.533) = 0.9 == 90%, a one-sided or one-tailed-test.



// http://faculty.vassar.edu/lowry/VassarStats.html
	// calculator gives mean = 37.8, difference from 'hypothetical population mean' == SRM 38.9 is -1.1, 
	// t = -1.975, 2 df, P 1-tailed = 0.093474, 2-tailed = 0.186948.

Example 1 using Student's t function.   ..\..\..\..\..\..\boost-sandbox\libs\mat
h_functions\example\students_t_example1.cpp Sat Jun 17 11:25:02 2006 140050727
mean = 37.8
standard deviation = 0.964365
 Student's t = -1.97566
Probability of t = -1.97566 is 0.0934283
Probability of t = 4.3 is 0.974971
Probability of t = 2.92 is 0.95
Probability of t = -2.92 is 0.0499996
Probability of Student's t = 0.05
Press any key to continue . . .




*/




