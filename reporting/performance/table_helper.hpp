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
#ifdef TEST_GSL
#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_version.h>

void gsl_handler(const char * reason, const char * file, int line, int gsl_errno)
{
   if(gsl_errno == GSL_ERANGE) return; // handle zero or infinity in our test code.
   throw std::domain_error(reason);
}

struct gsl_error_handler_setter
{
   gsl_error_handler_t * old_handler;
   gsl_error_handler_setter()
   {
      old_handler = gsl_set_error_handler(gsl_handler);
   }
   ~gsl_error_handler_setter()
   {
      gsl_set_error_handler(old_handler);
   }
};

static const gsl_error_handler_setter handler;

#endif

#ifdef TEST_RMATH
// Rmath overloads ftrunc, leading to strange errors from GCC unless we include this:
#include <boost/math/special_functions.hpp>
#define MATHLIB_STANDALONE
#include <Rmath.h>
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

