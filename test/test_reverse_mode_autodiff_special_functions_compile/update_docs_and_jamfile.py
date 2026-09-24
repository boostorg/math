#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#        Copyright Maksym Zhelyeznyakov 2025-2026
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
#           https://www.boost.org/LICENSE_1_0.txt)
#
import sys, os, subprocess
import logging
import re
from pathlib import Path


TEST_FOLDER="./"
SPECFUN_LIST=f"{TEST_FOLDER}specfun_list.txt"
CC="g++"
CPPFLAGS="--std=c++14"
CWD=os.getcwd()
BOOST_ROOT_DIR=f"{CWD}/../../../../"
BOOST_MATH_DIR=f"{BOOST_ROOT_DIR}/libs/math/"
BOOST_MATH_TEST_DIR=f"{BOOST_MATH_DIR}/test/"
INCLUDE_FLAGS=f"-I{BOOST_MATH_DIR}/include/ -I{BOOST_ROOT_DIR} -I{BOOST_MATH_TEST_DIR}"
COMPILATION_TABLE=f"{BOOST_ROOT_DIR}/libs/math/doc/differentiation/compilation_table.txt"
LOG_FILE=f"{CWD}/autogen.log"
LOG_FILE_VERBOSE=f"{CWD}/autogen_verbose.log"
JAM_COMPILE_COMMANDS_OUT="generated_jam_compile_commands.txt" 
DOC_FILE = f"{CWD}/../../doc/differentiation/autodiff_reverse.qbk"

def replace_between_markers(target_file, insert_file,
                            start_marker=r"\[/ BEGIN SPECFUN TABLE\]",
                            end_marker=r"\[/ END SPECFUN TABLE\]"):
    text = Path(target_file).read_text()
    replacement = Path(insert_file).read_text().rstrip("\n")
    pattern = re.compile(f"{start_marker}.*?{end_marker}", re.DOTALL)
    new_text = pattern.sub(f"{start_marker}\n{replacement}\n{end_marker}", text)
    Path(target_file).write_text(new_text)

if __name__=="__main__":
    replace_between_markers(DOC_FILE,COMPILATION_TABLE)
    replace_between_markers("../Jamfile.v2",JAM_COMPILE_COMMANDS_OUT,
                            start_marker="# BEGIN AUTODIFF LONG RUNNING TESTS",
                            end_marker="# END AUTODIFF LONG RUNNING TESTS")
