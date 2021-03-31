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
  sanitizer_test = [ "special_fun", "distribution_tests", "mp", "misc", "interpolators", "quadrature", "autodiff", "float128_tests" ]
  gnu_6_stds = [ "gnu++11", "gnu++14" ]
  gnu_8_stds = [ "gnu++11", "gnu++14", "gnu++17" ]
  gnu_10_stds = [ "gnu++11", "gnu++14", "gnu++17", "gnu++20" ]

  result = []

  for suite in things_to_test:
    #
    # Sanitizers:
    #
    result.append(linux_cxx("Ubunti g++-10 ASAN" + cxx + " " + suite, "g++-10", packages="g++-10", privileged=True, buildtype="boost", image="cppalliance/droneubuntu2004:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-10', 'CXXSTD': cxx, 'TEST_SUITE': suite, 'OPTIONS': '<cxxflags>-fsanitize=address <linkflags>-fsanitize=address <cxxflags>-DBOOST_CI_SANITIZER_BUILD' }, globalenv=globalenv))
    result.append(linux_cxx("Ubunti g++-10 USAN" + cxx + " " + suite, "g++-10", packages="g++-10", privileged=True, buildtype="boost", image="cppalliance/droneubuntu2004:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-10', 'CXXSTD': cxx, 'TEST_SUITE': suite, 'OPTIONS': '<cxxflags>-fsanitize=undefined <linkflags>-fsanitize=undefined' }, globalenv=globalenv))
    result.append(linux_cxx("Ubunti g++-10 TSAN" + cxx + " " + suite, "g++-10", packages="g++-10", privileged=True, buildtype="boost", image="cppalliance/droneubuntu2004:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-10', 'CXXSTD': cxx, 'TEST_SUITE': suite, 'OPTIONS': '<cxxflags>-fsanitize=thread <linkflags>-fsanitize=thread' }, globalenv=globalenv))


  return result

# from https://github.com/boostorg/boost-ci
load("@boost_ci//ci/drone/:functions.star", "linux_cxx","windows_cxx","osx_cxx","freebsd_cxx")
