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
  sanitizer_test = [ "special_fun", "distribution_tests", "misc", "interpolators", "quadrature", "float128_tests" ]
  gnu_5_stds = [ "gnu++14", "c++14" ]
  gnu_6_stds = [ "gnu++14", "c++14", "gnu++17", "c++17" ]
  clang_6_stds = [ "c++14", "c++17" ]
  gnu_9_stds = [ "gnu++17" ]
  clang_10_stds = [ "c++14", "c++17", "c++2a" ]

  result = []

  #for suite in sanitizer_test:
    #
    # Sanitizers:
    #
    #result.append(linux_cxx("Ubuntu g++-10 C++2a ASAN" + " " + suite, "g++-10", packages="g++-10", privileged=True, buildtype="boost", image="cppalliance/droneubuntu2004:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-10', 'CXXSTD': 'gnu++2a', 'TEST_SUITE': suite, 'OPTIONS': '<cxxflags>-fsanitize=address <linkflags>-fsanitize=address <cxxflags>-DBOOST_CI_SANITIZER_BUILD' }, globalenv=globalenv))
    #result.append(linux_cxx("Ubuntu g++-10 C++2a USAN" + " " + suite, "g++-10", packages="g++-10", privileged=True, buildtype="boost", image="cppalliance/droneubuntu2004:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-10', 'CXXSTD': 'gnu++2a', 'TEST_SUITE': suite, 'OPTIONS': '<cxxflags>-fsanitize=undefined <linkflags>-fsanitize=undefined <cxxflags>-DBOOST_CI_SANITIZER_BUILD' }, globalenv=globalenv))
    #result.append(linux_cxx("Ubuntu g++-10 C++2a TSAN" + " " + suite, "g++-10", packages="g++-10", privileged=True, buildtype="boost", image="cppalliance/droneubuntu2004:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-10', 'CXXSTD': 'gnu++2a', 'TEST_SUITE': suite, 'OPTIONS': '<cxxflags>-fsanitize=thread <linkflags>-fsanitize=thread <cxxflags>-DBOOST_CI_SANITIZER_BUILD' }, globalenv=globalenv))
    #result.append(linux_cxx("Ubuntu clang++-10 C++2a ISAN" + " " + suite, "clang++-10", packages="clang-10", privileged=True, buildtype="boost", image="cppalliance/droneubuntu2004:1", environment={'TOOLSET': 'clang', 'COMPILER': 'clang++-10', 'CXXSTD': 'gnu++2a', 'TEST_SUITE': suite, 'OPTIONS': '<cxxflags>-fsanitize=integer <linkflags>-fsanitize=integer' }, globalenv=globalenv))

  for suite in things_to_test:
    for cxx in gnu_9_stds:
      result.append(linux_cxx("Ubuntu g++ ARM64" + cxx + " " + suite, "g++", packages="g++", buildtype="boost", image="cppalliance/droneubuntu2204:multiarch", arch="arm64", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))

  return result

# from https://github.com/boostorg/boost-ci
load("@boost_ci//ci/drone/:functions.star", "linux_cxx","windows_cxx","osx_cxx","freebsd_cxx")
