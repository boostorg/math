// statistics_functions_example1.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Example 2 of using Thorsten Ottosen's statistics function.


//#include <boost/math/special_functions/students_t.hpp>
//	using boost::math::students_t;  // Probability of students_t(df, t).

//#define _SCL_SECURE_NO_DEPRECATE = 1 // avoid C4996 warning.
#define _SCL_SECURE_NO_DEPRECATE = 0 // get C4996 warning.

#ifdef _MSC_VER// needed for Boost.Array using const double
#  pragma warning(disable: 4510) // default constructor could not be generated.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
#endif

#include <iostream>
	using std::cout;
	using std::endl;
#include <iterator>
	using std::ostream_iterator;
#include <iomanip>
	using std::setprecision;
	using std::setw;
#include <cmath>
	using std::sqrt;
#include <vector>
	using std::vector;
#include <algorithm>
	using std::copy;
#include <list>
	using std::list;

#include <boost/array.hpp>
	using boost::array;
#include <boost/stat/univariate.hpp> // Various basic statistics functions.
	using boost::algorithm::mean;

#include <boost/math/special_functions/students_t.hpp>
	using boost::math::students_t;  // Probability of students_t(df, t).

int main()
{
	cout << "Example 1 of TO's statistics functions. ";
#if defined(__FILE__) && defined(__TIMESTAMP__)
	cout << "  " << __FILE__ << ' ' << __TIMESTAMP__ << ' '<< _MSC_FULL_VER;
#endif
	cout << endl;

	const int n = 5;
	double data [n] = {10.08, 10.11, 10.09, 10.10, 10.12}; // C array.
	//const array <double, n> adata = {10.08, 10.11, 10.09, 10.10, 10.12}; // Boost array c3892 connot assign to variable that is const
	array <const double, n> cddata = {10.08, 10.11, 10.09, 10.10, 10.12}; // Boost array C4510 default constructor could not be generated with const double
	//const array <const double, n> adata = {10.08, 10.11, 10.09, 10.10, 10.12}; // Boost array c4510 & c4512 assignment operator could not be generated
	// warning C4610: class 'boost::array<T,N>' can never be instantiated - user defined constructor required.
	// const array <const double, n> adata = {10.08, 10.11, 10.09, 10.10, 10.12}; // Boost array 
	//array <double, n> adata = {10.08, 10.11, 10.09, 10.10, 10.12}; // Boost array
	array <double, n> adata ; // Boost array
	copy (&data[0], &data[n], adata.begin()); // Note [n] (not [n-1]) because last is 1 beyond the end.
	// copy (&data[0], &data[n], adata.begin()); can't use this if double is const!
	// Error C3892: '_Dest' : you cannot assign to a variable that is const
	copy (cddata.begin(), cddata.end(), ostream_iterator<double>(cout, " "));	cout << endl;
	double cdm = mean<double>(cddata.begin(), cddata.end());
	BOOST_ASSERT(cddata.size() == n);
	cout <<"array size is " << cddata.size() << endl;
	cout << "array mean is " << cdm << endl;
	
	copy (adata.begin(), adata.end(), ostream_iterator<double>(cout, " "));	cout << endl;
	double am = mean<double>(adata.begin(), adata.end());
	BOOST_ASSERT(adata.size() == n);
	cout <<"array size is " << adata.size() << endl;
	cout << "array mean is " << am << endl;

	vector <double> vdata; // Std vector
	//vector <const double> vdata; // Std vector
	// vdata = {10.08, 10.11, 10.09, 10.10, 10.12}; is NOT allowed :-((
	vdata.reserve(5);
	vdata.assign(&data[0], &data[5]);
	cout << "vdata size = " << vdata.size() << endl;
	BOOST_ASSERT(vdata.size() == n);
	//copy (adata.begin(), adata.end(), vdata.begin()); // Asserts "vector iterator not dereferencable".
	copy (vdata.begin(), vdata.end(), ostream_iterator<double>(cout, " "));	cout << endl;
  double vm = mean<double>(vdata.begin(), vdata.end());
	cout << "vector mean is " << vm << endl;

	using boost::algorithm::variance;
	using boost::algorithm::std_deviation;

  double vv = variance<double>(vdata.begin(), vdata.end());
	cout << "vector variance is " << vv << endl;

  double vsd = std_deviation<double>(vdata.begin(), vdata.end());
	cout << "vector std_deviation is " << vsd << endl;

	using boost::algorithm::sorted_median;
	using boost::algorithm::unsorted_median;

 // double vsm = unsorted_median<double>(vdata.begin(), vdata.end());
	//cout << "vector sorted_median is " << vsm << endl;

	// Using contain list doesn't yet work for me - asked Thorsten Ottosen.
	//list<double> ldata;
	//ldata.assign(&data[0], &data[5]);
	//BOOST_ASSERT(ldata.size() == n);
	//copy (ldata.begin(), ldata.end(), ostream_iterator<double>(cout, " "));	cout << endl;
 // double lm = mean<double>(ldata.begin(), ldata.end());
	//cout << "list mean is " << lm << endl;

	double standard = 10.11;
	double t = (mean<double>(vdata.begin(), vdata.end()) - standard)
	* std::sqrt(static_cast<double>(n))
	/ std_deviation<double>(vdata.begin(), vdata.end());
	cout << "Student's t = " << t << endl; // Student's t = -1.41421
	double degrees_of_freedom = n-1;
	cout.precision(5); // Useful accuracy is only a few decimal digits, but seems to give at least 5.
	cout << "Probability of Student's t is " << students_t(degrees_of_freedom, abs(t)) << endl; // is 0.8849 , is 1 tailed.
	// So there is insufficient evidence of a difference to meet a 95% (1 in 20) criterion.



	return 0;
}  // int main()

/*

Output is:


*/
