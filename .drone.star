# Use, modification, and distribution are
# subject to the Boost Software License, Version 1.0. (See accompanying
# file LICENSE.txt)
#
# Copyright Rene Rivera 2020.
# Copyright John Maddock 2021.

# For Drone CI we use the Starlark scripting language to reduce duplication.
# As the yaml syntax for Drone CI is rather limited.
#
#
globalenv={}
linuxglobalimage="cppalliance/droneubuntu1604:1"
windowsglobalimage="cppalliance/dronevs2019"

def main(ctx):

  things_to_test = [ "special_fun", "distribution_tests", "mp", "misc", "interpolators", "quadrature", "autodiff", "long-running-tests", "float128_tests" ]
  nonx86_things_to_test = [ "special_fun", "distribution_tests", "mp", "misc", "interpolators", "quadrature", "autodiff", "long-running-tests" ]
  sanitizer_test = [ "special_fun", "distribution_tests", "misc", "interpolators", "quadrature", "float128_tests" ]
  gnu_5_stds = [ "gnu++11", "c++11", "gnu++14", "c++14" ]
  gnu_6_stds = [ "gnu++11", "c++11", "gnu++14", "c++14", "gnu++17", "c++17" ]
  clang_6_stds = [ "c++11", "c++14", "c++17" ]
  gnu_9_stds = [ "gnu++14", "c++14", "gnu++17", "c++17", "gnu++2a", "c++2a" ]
  clang_10_stds = [ "c++14", "c++17", "c++2a" ]
  nonx86_tests = [ "c++17" ]

  result = []

  for suite in nonx86_things_to_test:
    for cxx in nonx86_tests:
      result.append(osx_cxx("M1 Clang " + cxx + " " + suite, "clang++", buildscript="drone", buildtype="boost", xcode_version="14.1", environment={'TOOLSET': 'clang', 'CXXSTD': cxx, 'TEST_SUITE': suite, 'DEFINE': 'BOOST_MATH_NO_REAL_CONCEPT_TESTS,BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS,BOOST_MATH_MULTI_ARCH_CI_RUN', }, globalenv=globalenv))
      result.append(linux_cxx("ARM64: GCC 11" + cxx + " " + suite, "g++-11", packages="g++-11", buildscript="drone", buildtype="boost", image="cppalliance/droneubuntu2004:multiarch", environment={'TOOLSET': 'gcc', 'CXXSTD': cxx, 'TEST_SUITE': suite, 'DEFINE': 'BOOST_MATH_NO_REAL_CONCEPT_TESTS,BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS,BOOST_MATH_MULTI_ARCH_CI_RUN', }, arch="arm64", globalenv=globalenv))
      result.append(linux_cxx("S390X: GCC 11" + cxx + " " + suite, "g++-11", packages="g++-11", buildscript="drone", buildtype="boost", image="cppalliance/droneubuntu2004:multiarch", environment={'TOOLSET': 'gcc', 'CXXSTD': cxx, 'TEST_SUITE': suite, 'DEFINE': 'BOOST_MATH_NO_REAL_CONCEPT_TESTS,BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS,BOOST_MATH_MULTI_ARCH_CI_RUN', }, arch="s390x", globalenv=globalenv))


  return result

# from https://github.com/boostorg/boost-ci
load("@boost_ci//ci/drone/:functions.star", "linux_cxx","windows_cxx","osx_cxx","freebsd_cxx")
