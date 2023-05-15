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
  gcc13_things_to_test = [ "special_fun", "distribution_tests", "mp", "misc", "interpolators", "quadrature", "autodiff", "long-running-tests", "float128_tests", "new_floats" ]
  sanitizer_test = [ "special_fun", "distribution_tests", "misc", "interpolators", "quadrature", "float128_tests" ]
  gnu_5_stds = [ "gnu++14", "c++14" ]
  gnu_6_stds = [ "gnu++14", "c++14", "gnu++17", "c++17" ]
  clang_6_stds = [ "c++14", "c++17" ]
  gnu_9_stds = [ "gnu++14", "c++14", "gnu++17", "c++17", "gnu++2a", "c++2a" ]
  clang_10_stds = [ "c++14", "c++17", "c++2a" ]
  gnu_non_native = [ "gnu++17" ]
  gcc13_stds = [ "c++23" ]

  result = []

  for suite in things_to_test:
    for cxx in gcc13_stds:
      result.append(linux_cxx("Ubuntu g++-13 " + cxx + " " + suite, "g++-13", packages="g++-13", buildtype="boost", image="cppalliance/droneubuntu2304:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-13', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))

  return result

# from https://github.com/boostorg/boost-ci
load("@boost_ci//ci/drone/:functions.star", "linux_cxx","windows_cxx","osx_cxx","freebsd_cxx")
