# Use, modification, and distribution are
# subject to the Boost Software License, Version 1.0. (See accompanying
# file LICENSE.txt)
#
# Copyright Rene Rivera 2020.
# Copyright John Maddock 2021.

#!/bin/bash

set -ex
export TRAVIS_BUILD_DIR=$(pwd)
export DRONE_BUILD_DIR=$(pwd)
export TRAVIS_BRANCH=$DRONE_BRANCH
export VCS_COMMIT_ID=$DRONE_COMMIT
export GIT_COMMIT=$DRONE_COMMIT
export PATH=~/.local/bin:/usr/local/bin:$PATH

echo '==================================> BEFORE_INSTALL'

. .drone/before-install.sh

echo '==================================> INSTALL'

cd ..
if [ "$DRONE_BRANCH" == "master" ] || [[ "$DRONE_BRANCH" == */master ]]; then
    export BOOST_BRANCH="master"
else
    export BOOST_BRANCH="develop"
fi
git clone -b $BOOST_BRANCH --depth 1 https://github.com/boostorg/boost.git boost-root
cd boost-root
git submodule update --init tools/build
git submodule update --init libs/config
git submodule update --init libs/polygon
git submodule update --init tools/boost_install
git submodule update --init libs/headers
git submodule update --init tools/boostdep
cp -r $TRAVIS_BUILD_DIR/* libs/math
python tools/boostdep/depinst/depinst.py math
./bootstrap.sh
./b2 headers

if [[ $(uname) == "Linux" ]]; then
    echo 0 | sudo tee /proc/sys/kernel/randomize_va_space
fi

echo '==================================> BEFORE_SCRIPT'

. $DRONE_BUILD_DIR/.drone/before-script.sh

echo '==================================> SCRIPT'

# Require fused multiply-add, and let the compiler use it. gcc and clang both contract a*b + c by
# default, so -mfma is the only target flag needed. -O1 because gcc never fuses at -O0, and even at
# -O1 only with -fexpensive-optimizations, which enables the pass that forms FMAs. clang fuses
# within an expression at any level and rejects that flag.
case $TOOLSET in
    gcc*) FMA_FLAGS="<cxxflags>-O1 <cxxflags>-fexpensive-optimizations" ;;
    *) FMA_FLAGS="<cxxflags>-O1" ;;
esac
case $(uname -m) in
    x86_64)
        grep -qw fma /proc/cpuinfo || { echo "This runner's CPU has no FMA"; exit 1; }
        FMA_FLAGS="<cxxflags>-mfma $FMA_FLAGS"
        ;;
    *)
        # Fused multiply-add is part of the base aarch64 and s390x instruction sets: no -m flag needed.
        ;;
esac
echo "FMA: $(uname -m); $FMA_FLAGS"
echo "using $TOOLSET : : $COMPILER : <cxxflags>-std=$CXXSTD $OPTIONS $FMA_FLAGS ;" > ~/user-config.jam
(cd libs/config/test && ../../../b2 print_config_info print_math_info toolset=$TOOLSET)
(cd libs/math/test && ../../../b2 -d0 -j3 toolset=$TOOLSET $TEST_SUITE)

echo '==================================> AFTER_SUCCESS'

. $DRONE_BUILD_DIR/.drone/after-success.sh
