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

  things_to_test = [ "special_fun", "distribution_tests", "mp", "misc", "interpolators", "quadrature", "autodiff", "long-running-tests", "float128_tests", "concepts" ]
  gcc13_things_to_test = [ "special_fun", "distribution_tests", "mp", "misc", "interpolators", "quadrature", "autodiff", "long-running-tests", "float128_tests", "concepts", "new_floats" ]
  sanitizer_test = [ "special_fun", "distribution_tests", "misc", "interpolators", "quadrature", "float128_tests" ]
  reverse_mode_autodiff_test = [ "test_reverse_mode_autodiff", "autodiff-long-running-tests", "gradient_based_optimizers" ]
  cxx17_stds = [ "c++17" ]
  gnu_cxx17_stds = [ "gnu++17" ]
  non_native_stds = [ "c++17" ]
  gcc13_stds = [ "c++23" ]

  result = []

  for suite in sanitizer_test:
    #
    # Sanitizers:
    #
    result.append(linux_cxx("Ubuntu Clang-18 C++20 ASAN" + " " + suite, "clang++-18", packages="clang-18", privileged=True, buildtype="boost", image="cppalliance/droneubuntu2404:1", environment={'TOOLSET': 'clang', 'COMPILER': 'clang++-18', 'CXXSTD': 'gnu++20', 'TEST_SUITE': suite, 'OPTIONS': '<cxxflags>-fsanitize=address <linkflags>-fsanitize=address <cxxflags>-DBOOST_CI_SANITIZER_BUILD' }, globalenv=globalenv))
    result.append(linux_cxx("Ubuntu Clang-18 C++20 USAN" + " " + suite, "clang++-18", packages="clang-18", privileged=True, buildtype="boost", image="cppalliance/droneubuntu2404:1", environment={'TOOLSET': 'clang', 'COMPILER': 'clang++-18', 'CXXSTD': 'gnu++20', 'TEST_SUITE': suite, 'OPTIONS': '<cxxflags>-fsanitize=undefined <linkflags>-fsanitize=undefined <cxxflags>-DBOOST_CI_SANITIZER_BUILD' }, globalenv=globalenv))

  for suite in things_to_test:
    # Minimum supported C++17 compilers.
    for cxx in cxx17_stds:
      result.append(linux_cxx("Ubuntu g++-9 " + cxx + " " + suite, "g++-9", packages="g++-9", buildtype="boost", image="cppalliance/droneubuntu2004:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-9', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))
      result.append(linux_cxx("Ubuntu clang++-7 " + cxx + " " + suite, "clang++-7", packages="clang-7", llvm_os="bionic", llvm_ver="7", buildtype="boost", image="cppalliance/droneubuntu1804:1", environment={'TOOLSET': 'clang', 'COMPILER': 'clang++-7', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))

    # Keep one GNU-mode GCC configuration across the full suite.  This
    # exercises GNU extensions used by the __float128/libquadmath path.
    for cxx in gnu_cxx17_stds:
      result.append(linux_cxx("Ubuntu g++-12 " + cxx + " " + suite, "g++-12", packages="g++-12", buildtype="boost", image="cppalliance/droneubuntu2204:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-12', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))

    # Keep architecture coverage, but do it once in the minimum language mode.
    for cxx in non_native_stds:
      result.append(linux_cxx("Ubuntu g++ s390s " + cxx + " " + suite, "g++", packages="g++", buildtype="boost", image="cppalliance/droneubuntu2404:multiarch", arch="s390x", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))
      result.append(linux_cxx("Ubuntu g++ ARM64" + cxx + " " + suite, "g++", packages="g++", buildtype="boost", image="cppalliance/droneubuntu2404:multiarch", arch="arm64", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))
  for suite in gcc13_things_to_test:
    for cxx in gcc13_stds:
      result.append(linux_cxx("Ubuntu g++-14 " + cxx + " " + suite, "g++-14", packages="g++-14", buildtype="boost", image="cppalliance/droneubuntu2404:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-14', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))
  for suite in reverse_mode_autodiff_test:
    for cxx in gcc13_stds:
      result.append(linux_cxx("Ubuntu g++-14 " + cxx + " " + suite, "g++-14", packages="g++-14", buildtype="boost", image="cppalliance/droneubuntu2404:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-14', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))

  return result

# from https://github.com/boostorg/boost-ci
load("@boost_ci//ci/drone/:functions.star", "linux_cxx","windows_cxx","osx_cxx","freebsd_cxx")
