//  Copyright John Maddock 2007.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// Note this header must NOT include any other headers, for its
// use to be meaningful (because we use it in tests designed to
// detect missing includes).
//

static const float f = 0;
static const double d = 0;
static const long double l = 0;

template <class T>
void check_result_imp(T, T){}

template <class T1, class T2>
void check_result_imp(T1, T2)
{
   typedef int static_assertion[sizeof(T1) == 0xFFFF];
}

template <class T1, class T2>
void check_result(T2)
{
   T1 a = 0;
   T2 b = 0;
   return check_result_imp(a, b);
}
