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
  gnu_6_stds = [ "gnu++11", "gnu++14" ]
  gnu_8_stds = [ "gnu++11", "gnu++14", "gnu++17" ]
  gnu_10_stds = [ "gnu++11", "gnu++14", "gnu++17", "gnu++20" ]

  result = []

  for suite in things_to_test:
    for cxx in gnu_6_stds:
      result.append(linux_cxx("Ubunti g++-6 " + cxx + " " + suite, "g++-6", packages="g++-6", buildtype="boost", image=linuxglobalimage, environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-6', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))
    for cxx in gnu_8_stds:
      result.append(linux_cxx("Ubunti g++-8 " + cxx + " " + suite, "g++-8", packages="g++-8", buildtype="boost", image=linuxglobalimage, environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-8', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))
      result.append(linux_cxx("Ubunti g++-9 " + cxx + " " + suite, "g++-9", packages="g++-9", buildtype="boost", image=linuxglobalimage, environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-9', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))
      result.append(linux_cxx("Ubunti clang++-9 " + cxx + " " + suite, "clang++-9", packages="clang-9", llvm_os="xenial", llvm_ver="9", buildtype="boost", image=linuxglobalimage, environment={'TOOLSET': 'clang', 'COMPILER': 'clang++-9', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))
    for cxx in gnu_10_stds:
      result.append(linux_cxx("Ubunti g++-10 " + cxx + " " + suite, "g++-10", packages="g++-10", buildtype="boost", image="cppalliance/droneubuntu2004:1", environment={'TOOLSET': 'gcc', 'COMPILER': 'g++-10', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))
      result.append(linux_cxx("Ubunti clang++-10 " + cxx + " " + suite, "clang++-10", packages="clang-10", llvm_os="xenial", llvm_ver="10", buildtype="boost", image=linuxglobalimage, environment={'TOOLSET': 'clang', 'COMPILER': 'clang++-10', 'CXXSTD': cxx, 'TEST_SUITE': suite, }, globalenv=globalenv))

  return result

# from https://github.com/boostorg/boost-ci
load("@boost_ci//ci/drone/:functions.star", "linux_cxx","windows_cxx","osx_cxx","freebsd_cxx")
