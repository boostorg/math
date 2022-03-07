//  Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// The following change was released in VS 2022 17.1 Preview 1 and no longer allows compilation of *.c files
// https://github.com/microsoft/STL/pull/2148

#if !defined(_MSC_VER) || (_MSC_VER < 1930)
# error "MSVC1930 or newer run is NOT in effect".
#endif

int main()
{
   return 0;
}
