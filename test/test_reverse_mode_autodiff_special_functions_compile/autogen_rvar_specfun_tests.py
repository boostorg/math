#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# boost-no-inspect
#           Copyright Maksym Zhelyeznyakov 2025-2026
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
CPPFLAGS="--std=c++23"
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
logger = logging.getLogger("my_logger")
logger.setLevel(logging.INFO)
file_handler = logging.FileHandler(f"{LOG_FILE}", mode="w")
file_handler.setLevel(logging.INFO)
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.addHandler(console_handler)


logger2 = logging.getLogger("logger2")
logger2.setLevel(logging.INFO)
file_handler2 = logging.FileHandler(f"{LOG_FILE_VERBOSE}", mode="w")
file_handler2.setLevel(logging.INFO)
formatter2 = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler2.setFormatter(formatter2)
logger2.addHandler(file_handler2)
    
def generate_boost_test(func_sig,
                        et_str,
                        type_str,
                        min_val,
                        max_val,
                        cpp_file_dir):
    func_name = func_sig.split("(")[0]  # extract base function name
    
    filestr = f"""
    //           Copyright Maksym Zhelyeznyakov 2025-2026
    // Distributed under the Boost Software License, Version 1.0.
    // (See accompanying file LICENSE_1_0.txt or copy at
    //           https://www.boost.org/LICENSE_1_0.txt)
    //
    // This file was generated automatically with math/tests/autogen_rvar_specfun_tests.sh
    // DO NOT EDIT MANUALLY
    
    #define BOOST_MATH_REVERSE_MODE_ET_{et_str}
    #include <test_autodiff_reverse.hpp>
    #include <boost/math/special_functions.hpp>
    #include <boost/math/special_functions/logaddexp.hpp>
    #include <boost/math/tools/workaround.hpp>
    #include <cmath>
    
    BOOST_AUTO_TEST_SUITE(test_{func_name}_compiles)
    
    using namespace rdiff;
    
    BOOST_AUTO_TEST_CASE_TEMPLATE(test_{func_name}, T, {type_str})
    {{
        RandomSample<T> rng{{ {min_val}, {max_val} }};
        T               x = rng.next();
    
        rvar<T, 1> x_ad = x;
        auto y = boost::math::{func_sig.replace("arg", "x_ad")};
        auto y_expect = boost::math::{func_sig.replace("arg", "x")};
        BOOST_CHECK_CLOSE(y.item(), y_expect, 1000*boost_close_tol<T>());
    }}
    BOOST_AUTO_TEST_SUITE_END()
    """
    cpp_file_sig = f"test_reverse_mode_autodiff_{func_name}_compile_test_ET_{et_str}_{type_str}.cpp"
    cpp_file_loc = f"{cpp_file_dir}/{cpp_file_sig}"
    with open(cpp_file_loc,"w") as cpp_file:
        cpp_file.write(filestr)
    return cpp_file_loc

def compile_cpp_file(filename):
    base_file_name = filename.split("/")[-1]
    file_name_no_ext = base_file_name.split(".")[0]
    program_name = f"/tmp/{file_name_no_ext}"
    clean_up_command = f"rm {program_name}"
    compile_command = f"{CC} -o {program_name} {CPPFLAGS} {INCLUDE_FLAGS} -lboost_unit_test_framework {filename}"
    
    logger.info("cleaning up old executable")
    result_clean_up = subprocess.run(clean_up_command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    
    logger.info(f"attempting to compile {filename}")
    result_compile = subprocess.run(compile_command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    
    if result_compile.returncode == 0:
        logger.info(f"successfully compiled {filename}")
        return program_name, True
    else:
        logger.info(f"{filename} does not compile")
        logger2.error(f"{result_compile.stderr}")
        return program_name, False
    
if __name__=="__main__":
    bin_floats="bin_float_types"
    mp_type="multiprecision_float_types"
    
    floats_to_test = {
        "cpp_types" : bin_floats,
        "mp_types" : mp_type
        }

    compile_table_file = open(COMPILATION_TABLE, "w")
    jam_compile_commands = open(JAM_COMPILE_COMMANDS_OUT,"w")
    jam_compile_commands.write("test-suite autodiff-long-running-tests\n")
    jam_compile_commands.write("    :\n")
    compile_table_file.write("[table\n[[Function]\t[compiles with ET ON]\t[runs with ET ON]\t[compiles with ET OFF]	\t[runs with ET OFF]\t[works with multiprecision]\t[known issues]]\n")
    with open(SPECFUN_LIST,"r") as specfun_file:
        for line in specfun_file:
            if not line.count("#"):
                logger.info(line)
                split_str = line.split("\t")
                group, func_sig, minval, maxval = split_str[0], split_str[1], split_str[2], split_str[3]
                func_name = func_sig.split("(")[0] 
                if not os.path.exists(group):
                    os.makedirs(group)
                filename_et_on_cpp = generate_boost_test(func_sig, "ON", floats_to_test["cpp_types"], minval, maxval, group)
                filename_et_off_cpp = generate_boost_test(func_sig, "OFF", floats_to_test["cpp_types"], minval, maxval, group)
                filename_et_on_mp = generate_boost_test(func_sig, "ON", floats_to_test["mp_types"], minval, maxval, group)
                filename_et_off_mp = generate_boost_test(func_sig, "OFF", floats_to_test["mp_types"], minval, maxval, group)
        
                et_on_cpp_program, et_on_cpp_compiles = compile_cpp_file(filename_et_on_cpp)
                et_off_cpp_program, et_off_cpp_compiles = compile_cpp_file(filename_et_off_cpp)
                et_on_mp_program, et_on_mp_compiles = compile_cpp_file(filename_et_on_mp)
                et_off_mp_program, et_off_mp_compiles = compile_cpp_file(filename_et_off_mp)
                
                et_on_cpp_runs = "N/A"
                et_off_cpp_runs = "N/A"
                et_on_mp_runs = "N/A"
                et_off_mp_runs = "N/A"
                
                if et_on_cpp_compiles:
                    run_info = subprocess.run(f"{et_on_cpp_program}",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    if not run_info.returncode:
                        logger.info(f"{filename_et_on_cpp} ran successfully")
                        et_on_cpp_runs = "True"
                        jam_compile_commands.write(f"        [ run test_reverse_mode_autodiff_special_functions_compile/{filename_et_on_cpp} ]\n")
                    else:
                        logger.error(f"{filename_et_on_cpp} failed to run")
                        et_on_cpp_runs = "False"
                        jam_compile_commands.write(f"        [ run-fail test_reverse_mode_autodiff_special_functions_compile/{filename_et_on_cpp} ]\n")
                else:
                    jam_compile_commands.write(f"        [ compile-fail test_reverse_mode_autodiff_special_functions_compile/{filename_et_on_cpp} ]\n")
                    
                if et_off_cpp_compiles:
                    run_info = subprocess.run(f"{et_off_cpp_program}",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    if not run_info.returncode:
                        logger.info(f"{filename_et_off_cpp} ran successfully")
                        et_off_cpp_runs = "True"
                        jam_compile_commands.write(f"        [ run test_reverse_mode_autodiff_special_functions_compile/{filename_et_off_cpp} ]\n")
                    if not run_info.returncode:
                        logger.error(f"{filename_et_off_cpp} failed to run")
                        et_off_cpp_runs = "False"
                        jam_compile_commands.write(f"        [ run-fail test_reverse_mode_autodiff_special_functions_compile/{filename_et_off_cpp} ]\n")
                else:
                    jam_compile_commands.write(f"        [ compile-fail test_reverse_mode_autodiff_special_functions_compile/{filename_et_off_cpp} ]\n")
                
                if et_on_mp_compiles:
                    run_info = subprocess.run(f"{et_on_mp_program}",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    if not run_info.returncode:
                        logger.info(f"{filename_et_on_mp} ran successfully")
                        et_on_mp_runs = True
                        jam_compile_commands.write(f"        [ run test_reverse_mode_autodiff_special_functions_compile/{filename_et_on_mp} ]\n")
                    else:
                        logger.info(f"{filename_et_on_mp} failed to run")
                        et_on_mp_runs = False
                        jam_compile_commands.write(f"        [ run-fail test_reverse_mode_autodiff_special_functions_compile/{filename_et_on_mp} ]\n")
                else:
                    jam_compile_commands.write(f"        [ compile-fail test_reverse_mode_autodiff_special_functions_compile/{filename_et_on_mp} ]\n")
                    
                if et_off_mp_compiles:
                    run_info = subprocess.run(f"{et_off_mp_program}",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    if not run_info.returncode:
                        logger.info(f"{filename_et_off_mp} ran successfully")
                        et_off_mp_runs = True
                        jam_compile_commands.write(f"        [ run test_reverse_mode_autodiff_special_functions_compile/{filename_et_off_mp} ]\n")
                    if not run_info.returncode:
                        logger.error(f"{filename_et_off_mp} failed to run")
                        et_off_mp_runs = False
                        jam_compile_commands.write(f"        [ run-fail test_reverse_mode_autodiff_special_functions_compile/{filename_et_off_mp} ]\n")
                else:
                    jam_compile_commands.write(f"        [ compile-fail test_reverse_mode_autodiff_special_functions_compile/{filename_et_off_mp} ]\n")
                works_with_mp = et_on_mp_compiles and et_on_mp_runs and et_off_mp_compiles and et_off_mp_runs
                works_with_mp_str = ""
                if works_with_mp:
                    works_with_mp_str = "True"
                else:
                    works_with_mp = "False"
                    
                compile_table_file.write(f"[[{func_name}]\t[{et_on_cpp_compiles}]\t[{et_on_cpp_runs}]\t[{et_off_cpp_compiles}]\t[{et_off_cpp_runs}]\t[{works_with_mp}]\t[N/A]]\n")
    compile_table_file.write(
        """]"""
        )
    jam_compile_commands.write("    ;")
    compile_table_file.close()
    jam_compile_commands.close()
    
    
    
