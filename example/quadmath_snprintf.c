// Original source file copyright gcc.gnu.org, GNU Free Documentation License, Version 1.3.

// Notes for Boost.Math, (contains Quickbook snippets as C/C++ comments - do not remove!)
// Copyright Paul Bristow 2016.
// Distributed under the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


//[quadmath_snprintf_1
/*`Example of using GCC Quad-Precision Math Library quadmath `__float128` type,
taking a square root with sqrtq, and output using quadmath_snprintf.

From  GCC Quad-Precision Math Library,
[@https://gcc.gnu.org/onlinedocs/libquadmath.pdf 3.2 quadmath_snprintf, Convert to string]
(pages 9 and 10).

Requires GCC linker option `-lquadmath`.

If this linker option is missing then you will get errors like:
``
  /Cpp/float128/quadmath_snprintf/quadmath_snprintf.c:44: undefined reference to 'sqrtq'.
  /Cpp/float128/quadmath_snprintf/quadmath_snprintf.c:45: undefined reference to 'quadmath_snprintf'.
``
On one system, the header file (that contains all the `extern` declarations), included, for example:
``
  extern __float128 sqrtq (__float128) __quadmath_throw;
  extern __float128 strtoflt128 (const char *, char **) __quadmath_throw;
  extern int quadmath_snprintf (char *str, size_t size, const char *format, ...) __quadmath_throw;
``
An example of a location of `quadmath.h` is
``
  C:\program files\gcc-6-win64\lib\gcc\x86_64-w64-mingw32\6.1.1\include\quadmath.h
  ``
and library at
``
  C:\Program Files\gcc-6-win64\bin\libquadmath-0.dll
``

Command lines used (using [@http://www.codeblocks.org CodeBLocks]:
``
gcc.exe -Wall -g  -c J:\Cpp\float128\quadmath_snprintf\main.c -o obj\Debug\main.o
g++.exe  -o bin\Debug\quadmath_snprintf.exe obj\Debug\main.o  -lquadmath
``
*/
//] [/quadmath_snprintf_1]

#include <quadmath.h>
#include <stdlib.h>
#include <stdio.h>

int main ()
{
  __float128 r;
  int prec = 20;

  int width = 46;
  char buf[128];
  r = 2.0q;
  r = sqrtq (r);
  int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.20Qe", width, r);
  if ((size_t) n < sizeof buf)
  printf ("%s\n", buf);
  /* Prints: +1.41421356237309504880e+00 */
  quadmath_snprintf (buf, sizeof buf, "%Qa", r);
  if ((size_t) n < sizeof buf)
  printf ("%s\n", buf);
  /* Prints: 0x1.6a09e667f3bcc908b2fb1366ea96p+0 */
  n = quadmath_snprintf (NULL, 0, "%+-#46.*Qe", prec, r);
  if (n > -1)
  {
  char *str = malloc (n + 1);
  if (str)
  {
  quadmath_snprintf (str, n + 1, "%+-#46.*Qe", prec, r);
  printf ("%s\n", str);
  /* Prints: +1.41421356237309504880e+00 */
  }
  free (str);
  }
  return 0;
}

/*

Output:

+1.41421356237309504880e+00
0x1.6a09e667f3bcc908b2fb1366ea96p+0
+1.41421356237309504880e+00

*/

