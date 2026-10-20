#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# boost-no-inspect
#           Copyright Maksym Zhelyeznyakov 2025-2026
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
#           https://www.boost.org/LICENSE_1_0.txt)
#
import sys, os, subprocess, shutil, tempfile
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
    
    namespace test_{func_name}_compiles {{
    
    using namespace rdiff;
    
    template <typename T>
    void test_{func_name}()
    {{
        RandomSample<T> rng{{ {min_val}, {max_val} }};
        T               x = rng.next();
    
        rvar<T, 1> x_ad = x;
        auto y = boost::math::{func_sig.replace("arg", "x_ad")};
        auto y_expect = boost::math::{func_sig.replace("arg", "x")};
        CHECK_RELATIVE_ERROR(y.item(), y_expect, 1000*boost_close_tol<T>() / 100);
    }}
    }} // namespace test_{func_name}_compiles
    
    int main()
    {{
        mp11::mp_for_each<mp11::mp_transform<mp11::mp_identity, {type_str}>>([](auto t) {{ test_{func_name}_compiles::test_{func_name}<typename decltype(t)::type>(); }});
        return boost::math::test::report_errors();
    }}
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
    compile_command = f"{CC} -o {program_name} {CPPFLAGS} {INCLUDE_FLAGS} {filename}"
    
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
    
def check_and_keep(filename, group, jam_compile_commands):
    """Compile and run a generated test.  If it passes, copy it into `group`
    and list it in the Jamfile; otherwise make sure no stale copy is left in
    the repository.  Returns (compiles, runs) as table strings ("N/A" if the
    test did not compile)."""
    base_file_name = os.path.basename(filename)
    repo_file = f"{group}/{base_file_name}"
    program_name, compiles = compile_cpp_file(filename)
    runs = "N/A"
    if compiles:
        run_info = subprocess.run(f"{program_name}",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if not run_info.returncode:
            logger.info(f"{filename} ran successfully")
            runs = "True"
        else:
            logger.error(f"{filename} failed to run")
            runs = "False"
    if runs == "True":
        if not os.path.exists(group):
            os.makedirs(group)
        shutil.copyfile(filename, repo_file)
        jam_compile_commands.write(f"        [ run test_reverse_mode_autodiff_special_functions_compile/{repo_file} ]\n")
    elif os.path.exists(repo_file):
        os.remove(repo_file)
    return str(compiles), runs

if __name__=="__main__":
    bin_floats="bin_float_types"
    mp_type="multiprecision_float_types"
    
    floats_to_test = {
        "cpp_types" : bin_floats,
        "mp_types" : mp_type
        }

    scratch_dir = tempfile.mkdtemp()
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
                # Generate into a scratch directory; only tests that compile and
                # pass are copied into the repository and listed in the Jamfile.
                filename_et_on_cpp = generate_boost_test(func_sig, "ON", floats_to_test["cpp_types"], minval, maxval, scratch_dir)
                filename_et_off_cpp = generate_boost_test(func_sig, "OFF", floats_to_test["cpp_types"], minval, maxval, scratch_dir)
                filename_et_on_mp = generate_boost_test(func_sig, "ON", floats_to_test["mp_types"], minval, maxval, scratch_dir)
                filename_et_off_mp = generate_boost_test(func_sig, "OFF", floats_to_test["mp_types"], minval, maxval, scratch_dir)

                et_on_cpp_compiles, et_on_cpp_runs = check_and_keep(filename_et_on_cpp, group, jam_compile_commands)
                et_off_cpp_compiles, et_off_cpp_runs = check_and_keep(filename_et_off_cpp, group, jam_compile_commands)
                et_on_mp_compiles, et_on_mp_runs = check_and_keep(filename_et_on_mp, group, jam_compile_commands)
                et_off_mp_compiles, et_off_mp_runs = check_and_keep(filename_et_off_mp, group, jam_compile_commands)

                works_with_mp = et_on_mp_runs == "True" and et_off_mp_runs == "True"
                works_with_mp_str = "True" if works_with_mp else "False"

                compile_table_file.write(f"[[{func_name}]\t[{et_on_cpp_compiles}]\t[{et_on_cpp_runs}]\t[{et_off_cpp_compiles}]\t[{et_off_cpp_runs}]\t[{works_with_mp_str}]\t[N/A]]\n")
    compile_table_file.write(
        """]"""
        )
    jam_compile_commands.write("    ;")
    compile_table_file.close()
    jam_compile_commands.close()
    shutil.rmtree(scratch_dir)
    
    
    
