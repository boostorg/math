//  Copyright John Maddock 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef TABLE_HELPER_HPP
#define TABLE_HELPER_HPP

#include <vector>
#include <string>

//
// Also include headers for whatever else we may be testing:
//
#ifdef TEST_LIBSTDCXX
#include <tr1/cmath>
#include <stdexcept>
#endif

extern std::vector<std::vector<double> > data;

std::string sanitize_string(const std::string& s);
void load_table(std::vector<std::vector<std::string> >& table, std::string::const_iterator begin, std::string::const_iterator end);
std::string save_table(std::vector<std::vector<std::string> >& table);
void add_to_all_sections(const std::string& id, std::string list_name = "all_sections");
void add_cell(const std::string& cell_name, const std::string& table_name, const std::string& row_name, const std::string& column_heading);
void set_result(const std::string& cell_name, const std::string& cell_content, const std::string& table_name, const std::string& row_name, const std::string& column_name);
void report_execution_time(double t, std::string table, std::string row, std::string heading);
std::string get_compiler_options_name();

#endif // TABLE_HELPER_HPP

