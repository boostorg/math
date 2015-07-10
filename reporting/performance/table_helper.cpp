//  Copyright John Maddock 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning (disable : 4224)
#endif

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/interprocess/sync/named_mutex.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <vector>
#include <set>
#include <iostream>
#include "table_helper.hpp"

std::vector<std::vector<double> > data;

inline std::string sanitize_string(const std::string& s)
{
   static const boost::regex e("[^a-zA-Z0-9]+");
   std::string result = boost::regex_replace(s, e, "_");
   while(result[0] == '_')
      result.erase(0);
   return result;
}

static std::string content;
boost::filesystem::path path_to_content;

struct content_loader
{
   boost::interprocess::named_mutex mu;
   boost::interprocess::scoped_lock<boost::interprocess::named_mutex> lock;
   content_loader() : mu(boost::interprocess::open_or_create, "handle_test_result"), lock(mu)
   {
      boost::filesystem::path p(__FILE__);
      p = p.parent_path();
      p /= "doc";
      p /= "performance_tables.qbk";
      path_to_content = p;
      if(boost::filesystem::exists(p))
      {
         boost::filesystem::ifstream is(p);
         if(is.good())
         {
            do
            {
               char c = static_cast<char>(is.get());
               if(c != EOF)
                  content.append(1, c);
            } while(is.good());
         }
      }
   }
   ~content_loader()
   {
      boost::filesystem::ofstream os(path_to_content);
      os << content;
   }
   void instantiate()const
   {
   }
};

static const content_loader loader;

void load_table(std::vector<std::vector<std::string> >& table, std::string::const_iterator begin, std::string::const_iterator end)
{
   static const boost::regex item_e(
      "\\["
      "([^\\[\\]]*(?0)?)*"
      "\\]"
      );

   boost::regex_token_iterator<std::string::const_iterator> i(begin, end, item_e), j;

   while(i != j)
   {
      // Add a row:
      table.push_back(std::vector<std::string>());
      boost::regex_token_iterator<std::string::const_iterator> k(i->first + 1, i->second - 1, item_e);
      while(k != j)
      {
         // Add a cell:
         table.back().push_back(std::string(k->first + 1, k->second - 1));
         ++k;
      }
      ++i;
   }
}

std::string save_table(std::vector<std::vector<std::string> >& table)
{
   std::string result;

   static const boost::regex value_e("\\d+");

   std::vector<boost::uintmax_t> best_values(table[0].size(), (std::numeric_limits<boost::uintmax_t>::max)());

   for(unsigned i = 1; i < table[0].size(); ++i)
   {
      for(unsigned k = 1; k < table.size(); ++k)
      {
         boost::smatch what;
         if(regex_search(table[k][i], what, value_e))
         {
            boost::uintmax_t val = boost::lexical_cast<boost::uintmax_t>(what.str());
            if(val < best_values[i])
               best_values[i] = val;
         }
      }
   }

   for(std::vector<std::vector<std::string> >::const_iterator i = table.begin(), j = table.end(); i != j; ++i)
   {
      result += "[";
      for(std::vector<std::string>::const_iterator k = i->begin(), l = i->end(); k != l; ++k)
      {
         result += "[";
         if((i == table.begin()) || (k == i->begin()))
         {
            result += *k;
         }
         else
         {
            boost::smatch what;
            if(regex_search(*k, what, value_e))
            {
               boost::uintmax_t val = boost::lexical_cast<boost::uintmax_t>(what.str());
               bool have_role = false;
               if(val < 1.2 * best_values[k - i->begin()])
               {
                  have_role = true;
                  result += "[role green ";
               }
               else if(val > 2 * best_values[k - i->begin()])
               {
                  have_role = true;
                  result += "[role red ";
               }
               result += boost::lexical_cast<std::string>(val);
               result += "ns";
               if(have_role)
                  result += "]";
            }
            else
               result += *k;
         }
         result += "]";
      }
      result += "]\n";
   }
   return result;
}

void add_to_all_sections(const std::string& id, std::string list_name)
{
   std::string::size_type pos = content.find("[template " + list_name + "[]"), end_pos;
   if(pos == std::string::npos)
   {
      //
      // Just append to the end:
      //
      content.append("\n[template ").append(list_name).append("[]\n[").append(id).append("]\n]\n");
   }
   else
   {
      //
      // Read in the all list of sections, add our new one (in alphabetical order),
      // and then rewrite the whole thing:
      //
      static const boost::regex item_e(
         "\\["
         "([^\\[\\]]*(?0)?)*"
         "\\]|\\]"
         );
      boost::regex_token_iterator<std::string::const_iterator> i(content.begin() + pos + 12 + list_name.size(), content.end(), item_e), j;
      std::set<std::string> sections;
      while(i != j)
      {
         if(i->length() == 1)
         {
            end_pos = i->first - content.begin();
            break;
         }
         sections.insert(std::string(i->first + 1, i->second - 1));
         ++i;
      }
      sections.insert(id);
      std::string new_list = "\n";
      for(std::set<std::string>::const_iterator sec = sections.begin(); sec != sections.end(); ++sec)
      {
         new_list += "[" + *sec + "]\n";
      }
      content.replace(pos + 12 + list_name.size(), end_pos - pos - 12 - list_name.size(), new_list);
   }
}

void add_cell(const std::string& cell_name, const std::string& table_name, const std::string& row_name, const std::string& column_heading)
{
   //
   // Load the table, add our data, and re-write:
   //
   std::string table_id = "table_" + sanitize_string(table_name);
   boost::regex table_e("\\[table:" + table_id
      + "\\s[^\\[]+"
      "((\\["
      "([^\\[\\]]*(?2)?)*"
      "\\]\\s*)*\\s*)"
      "\\]"
      );

   boost::smatch table_location;
   if(regex_search(content, table_location, table_e))
   {
      std::vector<std::vector<std::string> > table_data;
      load_table(table_data, table_location[1].first, table_location[1].second);
      //
      // Figure out which column we're on:
      //
      unsigned column_id = 1001u;
      for(unsigned i = 0; i < table_data[0].size(); ++i)
      {
         if(table_data[0][i] == column_heading)
         {
            column_id = i;
            break;
         }
      }
      if(column_id > 1000)
      {
         //
         // Need a new column, must be adding a new compiler to the table!
         //
         table_data[0].push_back(column_heading);
         for(unsigned i = 1; i < table_data.size(); ++i)
            table_data[i].push_back(std::string());
         column_id = table_data[0].size() - 1;
      }
      //
      // Figure out the row:
      //
      unsigned row_id = 1001;
      for(unsigned i = 1; i < table_data.size(); ++i)
      {
         if(table_data[i][0] == row_name)
         {
            row_id = i;
            break;
         }
      }
      if(row_id > 1000)
      {
         //
         // Need a new row, add it now:
         //
         table_data.push_back(std::vector<std::string>());
         table_data.back().push_back(row_name);
         for(unsigned i = 1; i < table_data[0].size(); ++i)
            table_data.back().push_back(std::string());
         row_id = table_data.size() - 1;
      }
      //
      // Update the entry:
      //
      std::string& s = table_data[row_id][column_id];
      if(s.empty())
      {
         std::cout << "Adding " << cell_name << " to empty cell.";
         s = "[" + cell_name + "]";
      }
      else
      {
         if(cell_name.find("_boost_") != std::string::npos)
         {
            std::cout << "Adding " << cell_name << " to start of cell.";
            s.insert(0, "[" + cell_name + "][br][br]");
         }
         else
         {
            std::cout << "Adding " << cell_name << " to end of cell.";
            if((s.find("_boost_") != std::string::npos) && (s.find("[br]") == std::string::npos))
               s += "[br]"; // extra break if we're adding directly after the boost results.
            s += "[br][" + cell_name + "]";
         }
      }
      //
      // Convert back to a string and insert into content:
      std::string c = save_table(table_data);
      content.replace(table_location.position(1), table_location.length(1), c);
   }
   else
   {
      //
      // Create a new table and try again:
      //
      std::string new_table = "\n[template " + table_id;
      new_table += "[]\n[table:" + table_id;
      new_table += " Performance comparison for ";
      new_table += table_name;
      new_table += "\n[[Function][";
      new_table += column_heading;
      new_table += "]]\n";
      new_table += "[[";
      new_table += row_name;
      new_table += "][[";
      new_table += cell_name;
      new_table += "]]]\n]\n]\n";

      std::string::size_type pos = content.find("[/tables:]");
      if(pos != std::string::npos)
         content.insert(pos + 10, new_table);
      else
         content += "\n\n[/tables:]\n" + new_table;
      //
      // Add a section for this table as well:
      //
      std::string section_id = "section_" + sanitize_string(table_name);
      if(content.find(section_id + "[]") == std::string::npos)
      {
         std::string new_section = "\n[template " + section_id + "[]\n[section:" + section_id + " " + table_name + "]\n[" + table_id + "]\n[endsect]\n]\n";
         pos = content.find("[/sections:]");
         if(pos != std::string::npos)
            content.insert(pos + 12, new_section);
         else
            content += "\n\n[/sections:]\n" + new_section;
         add_to_all_sections(section_id);
      }
      //
      // Add to list of all tables (not in sections):
      //
      add_to_all_sections(table_id, "all_tables");
   }
}

void set_result(const std::string& cell_name, const std::string& cell_content, const std::string& table_name, const std::string& row_name, const std::string& column_name)
{
   loader.instantiate();
   const boost::regex e("\\[template\\s+" + cell_name +
      "\\[\\]([^\\n]*)\\]$");

   boost::smatch what;
   if(regex_search(content, what, e))
   {
      content.replace(what.position(1), what.length(1), cell_content);
   }
   else
   {
      // Need to add new content:
      std::string::size_type pos = content.find("[/Cell Content:]");
      std::string t = "\n[template " + cell_name + "[] " + cell_content + "]";
      if(pos != std::string::npos)
         content.insert(pos + 16, t);
      else
      {
         content.insert(0, t);
         content.insert(0, "[/Cell Content:]");
      }
   }
   //
   // Check to verify that our content is actually used somewhere,
   // if not we need to create a place for it:
   //
   if(content.find("[" + cell_name + "]") == std::string::npos)
      add_cell(cell_name, table_name, row_name, column_name);
}


void report_execution_time_multi_compilation(double t, std::string function_group, std::string function, std::string compilation)
{
   std::string cell_name = sanitize_string(function) + sanitize_string(compilation) + sanitize_string(BOOST_COMPILER) + sanitize_string(BOOST_PLATFORM);
   std::string cell_content = boost::lexical_cast<std::string>(static_cast<boost::uintmax_t>(t / 1e-9));
   cell_content += "ns";
   std::string table_name = function_group;
   std::string row_name = function;
   std::string column_name = compilation;
   set_result(cell_name, cell_content, table_name, row_name, column_name);
}

std::string get_compiler_options_name()
{
#if defined(BOOST_MSVC)
   std::string result;
#ifdef _DEBUG
   result =  "cl /Od";
#elif defined(__AVX2__)
   result = "cl /arch:AVX2 /Ox";
#elif defined(__AVX__)
   result = "cl /arch:AVX /Ox";
#elif _M_IX86_FP == 2
   result = "cl /arch:sse2 /Ox";
#else
   result = "cl /arch:ia32 /Ox";
#endif
#ifdef _M_AMD64
   result += " (x64 build)";
#else
   result += " (x86 build)";
#endif
   return result;
#else
   return "Unknown";
#endif
}

