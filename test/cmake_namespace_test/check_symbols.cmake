# Copyright 2026 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# Usage: cmake -DNM=<nm> -DFILE=<binary> -P check_symbols.cmake
# Fails if any symbol is mangled into ::boost, or if nothing was emitted into ::my_lib::math.

execute_process(COMMAND "${NM}" "${FILE}" OUTPUT_VARIABLE symbols RESULT_VARIABLE result)

if(NOT result EQUAL 0)
  message(FATAL_ERROR "${NM} failed on ${FILE}")
endif()

string(REGEX MATCHALL "[^\n]*N5boost[^\n]*" leaked "${symbols}")
if(leaked)
  list(LENGTH leaked count)
  string(REPLACE ";" "\n" leaked "${leaked}")
  message(FATAL_ERROR "${count} symbols are still in namespace boost:\n${leaked}")
endif()

if(NOT symbols MATCHES "N6my_lib4math")
  message(FATAL_ERROR "No symbols found in namespace my_lib::math, the check is not meaningful")
endif()

message(STATUS "No symbols in namespace boost")
