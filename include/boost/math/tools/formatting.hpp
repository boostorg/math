//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_FORMATTING_HPP
#define BOOST_MATH_TOOLS_FORMATTING_HPP

#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <cmath>
#include <boost/throw_exception.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_detected.hpp>
#include <boost/type_traits/is_class.hpp>
#include <boost/type_traits/is_arithmetic.hpp>

#ifdef _MSC_VER
#pragma once
#endif

namespace boost {
   namespace math {
      namespace tools {

         template <class T, bool = boost::is_arithmetic<T>::value || boost::is_class<T>::value>
         struct is_integer_like
         {
            static const bool value = false;
         };
         template <class T>
         struct is_integer_like<T, true>
         {
            static const bool value = std::numeric_limits<T>::is_integer;
         };

         template <class T, bool = boost::is_arithmetic<T>::value || boost::is_class<T>::value>
         struct is_float_like
         {
            // The point "floats":
            static const bool value = false;
         };
         template <class T>
         struct is_float_like<T, true>
         {
            // The point "floats":
            static const bool value = std::numeric_limits<T>::max_exponent != std::numeric_limits<T>::min_exponent;
         };

         template<class T>
         using real_t = decltype(std::declval<T&>().real(std::declval<T const&>().imag()));
         template<class T>
         using imag_t = decltype(std::declval<T&>().imag(std::declval<T const&>().real()));
         template <class T>
         using polynomial_t = decltype(std::declval<const T&>()(std::declval<const T&>()[0] + std::declval<const T&>().order()));

         template <class T>
         struct is_complex_like
         {
            static const bool value = boost::is_detected_v<real_t, T> && boost::is_detected_v<imag_t, T>;
         };

         template <class T>
         struct is_unhandled_type
         {
            static const bool value = !is_complex_like<T>::value && !is_float_like<T>::value && !is_integer_like<T>::value;
         };


         enum output_format_t
         {
            text_format,
            docbook_format,
            latex_format,
            html_format
         };

         enum styling_level_t
         {
            no_styling = 0,
            minimal_styling = 1,
            full_styling = 6
         };

         enum multiplyer_t
         {
            multiply_times = 0,
            multiply_dot = 1,
            multiply_x = 2
         };

         enum imaginary_i_t
         {
            upright_i = 0,
            slanted_i = 1,
            doublestruck_i = 2
         };

         enum zero_component_t
         {
            hide_zero_components = 0,
            show_zero_components = 1
         };

         struct latex_as_text_t {};
         struct latex_as_equation_t {};

         constexpr const latex_as_text_t* latex_as_text = nullptr;
         constexpr const latex_as_equation_t* latex_as_equation = nullptr;

         template <output_format_t Format, class charT, class Traits>
         class basic_numeric_formatter;

         template <class charT, class Traits>
         class basic_numeric_formatter_base
         {
            styling_level_t styling_level;
            imaginary_i_t i_style;
            multiplyer_t  multiply_style;
            bool requires_parenthesis;
            bool m_show_zero_components;
         protected:
            template <class Float>
            static void decompose_float(std::basic_ostream<charT, Traits>& os, const Float& f, std::basic_string<charT, Traits>& mantissa, std::basic_string<charT, Traits>& exponent)
            {
               static const charT exponent_string[2] = { 'e', 'E' };
               std::basic_stringstream<charT, Traits> ss;
               ss.copyfmt(os);
               ss.imbue(os.getloc());
               ss << f;
               std::basic_string<charT, Traits> s = ss.str();
               typename std::basic_string<charT, Traits>::size_type pos = s.find_first_of(exponent_string);
               if (pos != std::basic_string<charT, Traits>::npos)
               {
                  mantissa.assign(s, 0, pos);
                  exponent.assign(s, pos + 1, s.size() - pos - 1);
               }
               else
                  mantissa = s;
            }
         public:
            basic_numeric_formatter_base() : styling_level(full_styling), requires_parenthesis(false), multiply_style(multiply_times), m_show_zero_components(false) {}

            void styling(styling_level_t i)
            {
               styling_level = i;
            }
            styling_level_t styling()const
            {
               return styling_level;
            }
            void imaginary_style(imaginary_i_t i)
            {
               i_style = i;
            }
            imaginary_i_t imaginary_style()const
            {
               return i_style;
            }
            void multiply(multiplyer_t t)
            {
               multiply_style = t;
            }
            multiplyer_t multiply()const
            {
               return multiply_style;
            }
            void parenthesis(bool b)
            {
               requires_parenthesis = b;
            }
            bool parenthesis()const
            {
               return requires_parenthesis;
            }
            void show_zero_components(bool b)
            {
               m_show_zero_components = b;
            }
            bool show_zero_components()const
            {
               return m_show_zero_components;
            }
            template <class Integer>
            static std::basic_ostream<charT, Traits>& format_integer(std::basic_ostream<charT, Traits>& os, const Integer& i)
            {
               return os << i;
            }
            template <class Float>
            static std::basic_ostream<charT, Traits>& format_float(std::basic_ostream<charT, Traits>& os, const Float& f)
            {
               std::basic_string<charT, Traits> mantissa, exponent;
               decompose_float(os, f, mantissa, exponent);
               if (exponent.size())
                  return os << mantissa << "x10^" << exponent;
               return os << mantissa;
            }
            template <class Complex>
            std::basic_ostream<charT, Traits>& format_complex(std::basic_ostream<charT, Traits>& os, const Complex& f)
            {
               if (!this->show_zero_components())
               {
                  if (f.imag() == 0)
                  {
                     format_float(os, f.real());
                     return os;
                  }
                  else if (f.real() == 0)
                  {
                     format_float(os, f.imag());
                     os.put(os.widen('i'));
                     return os;
                  }
               }
               format_float(os, f.real());
               typename Complex::value_type i(f.imag());
               bool isneg = i < 0;
               if (isneg)
               {
                  i = -i;
                  os << " - ";
               }
               else
                  os << " + ";
               format_float(os, i);
               os.put(os.widen('i'));
               return os;
            }
            bool latex_as_equation()const { return false; }
            void latex_as_equation(bool) {}
         };

         template <class charT, class Traits>
         class basic_numeric_formatter<text_format, charT, Traits> : public basic_numeric_formatter_base<charT, Traits>
         {
         };

         template <class charT, class Traits>
         class basic_numeric_formatter<docbook_format, charT, Traits> : public basic_numeric_formatter_base<charT, Traits>
         {
         protected:
            template <class Integer>
            std::basic_ostream<charT, Traits>& format_integer(std::basic_ostream<charT, Traits>& os, const Integer& i)
            {
               if (this->styling() >= minimal_styling)
                  os << "<phrase role=\"number\">";
               if (this->styling() > minimal_styling)
                  os << "<phrase role=\"integer\">";
               os << i;
               if (this->styling() > minimal_styling)
                  os << "</phrase>";
               if (this->styling() >= minimal_styling)
                  os << "</phrase>";
               return os;
            }
            template <class Float>
            std::basic_ostream<charT, Traits>& format_float(std::basic_ostream<charT, Traits>& os, const Float& f)
            {
               if (this->styling() >= minimal_styling)
                  os << "<phrase role=\"number\">";
               if (this->styling() > minimal_styling)
                  os << "<phrase role=\"float\">";

               if ((std::isinf)(f))
               {
                  if (f < 0)
                     os.put(os.widen('-'));
                  os << "&#8734;";
               }
               else if ((std::isnan(f)))
               {
                  os << "NaN";
               }
               else
               {
                  std::basic_string<charT, Traits> mantissa, exponent;
                  basic_numeric_formatter_base<charT, Traits>::decompose_float(os, f, mantissa, exponent);

                  if (exponent.size())
                  {
                     os << mantissa;
                     if (this->multiply() == multiply_x)
                        os.put(os.widen('x'));
                     else if (this->multiply() == multiply_dot)
                        os << "&#x22C5;";
                     else
                        os << "&#x00D7;";
                     os << "10<superscript>" << exponent << "</superscript>";
                  }
                  else
                     os << mantissa;
               }

               if (this->styling() > minimal_styling)
                  os << "</phrase>";
               if (this->styling() >= minimal_styling)
                  os << "</phrase>";

               return os;
            }
            template <class Complex>
            std::basic_ostream<charT, Traits>& format_complex(std::basic_ostream<charT, Traits>& os, const Complex& f)
            {
               if (this->styling() >= minimal_styling)
                  os << "<phrase role=\"number\">";
               if (this->styling() > minimal_styling)
                  os << "<phrase role=\"complex\">";

               if ((std::isnan(f.real())) || (std::isnan)(f.imag()))
               {
                  os << "NaN";
               }
               else if ((std::isinf)(f.real()) || (std::isinf)(f.imag()))
               {
                  os << "&#8734;&#x0303;";
               }
               else
               {
                  styling_level_t saved_style = this->styling();
                  this->styling(no_styling);

                  bool need_i = true;

                  if (!this->show_zero_components() && (f.imag() == 0))
                  {
                     format_float(os, f.real());
                     need_i = false;
                  }
                  else if (!this->show_zero_components() && (f.real() == 0))
                  {
                     format_float(os, f.imag());
                  }
                  else
                  {
                     format_float(os, f.real());
                     typename Complex::value_type i(f.imag());
                     bool isneg = i < 0;
                     if (isneg)
                     {
                        i = -i;
                        os << " - ";
                     }
                     else
                        os << " + ";
                     format_float(os, i);
                  }
                  if (need_i)
                  {
                     if (saved_style >= minimal_styling)
                        os << "<phrase role=\"imaginary_i\">";
                     if (this->imaginary_style() == doublestruck_i)
                        os << "&#x2148;";
                     else if (this->imaginary_style() == slanted_i)
                        os << "<emphasis>i</emphasis>";
                     else
                        os.put(os.widen('i'));
                     if (saved_style >= minimal_styling)
                        os << "</phrase>";
                  }
                  this->styling(saved_style);
               }

               if (this->styling() > minimal_styling)
                  os << "</phrase>";
               if (this->styling() >= minimal_styling)
                  os << "</phrase>";

               return os;
            }
         };

         template <class charT, class Traits>
         class basic_numeric_formatter<latex_format, charT, Traits> : public basic_numeric_formatter_base<charT, Traits>
         {
            bool is_latex_as_equation;
            bool inside_equation;
         public:
            basic_numeric_formatter() : is_latex_as_equation(true), inside_equation(false) {}
            void latex_as_equation(bool b) { is_latex_as_equation = b; }
            bool latex_as_equation()const { return is_latex_as_equation; }

            template <class Integer>
            std::basic_ostream<charT, Traits>& format_integer(std::basic_ostream<charT, Traits>& os, const Integer& i)
            {
               if (!inside_equation && latex_as_equation())
                  os.put(os.widen('$'));
               basic_numeric_formatter_base<charT, Traits>::format_integer(os, i);
               if (!inside_equation && latex_as_equation())
                  os.put(os.widen('$'));
               return os;
            }

            template <class Float>
            std::basic_ostream<charT, Traits>& format_float(std::basic_ostream<charT, Traits>& os, const Float& f)
            {
               std::basic_string<charT, Traits> mantissa, exponent;
               basic_numeric_formatter_base<charT, Traits>::decompose_float(os, f, mantissa, exponent);

               if ((std::isinf)(f))
               {
                  if (!inside_equation)
                     os.put(os.widen('$'));
                  if (f < 0)
                     os.put(os.widen('-'));
                  os << "\\infty";
                  if (!inside_equation)
                     os.put(os.widen('$'));
                  return os;
               }
               else if ((std::isnan)(f))
               {
                  if (inside_equation)
                     os << "\\mathrm{NaN}";
                  else
                     os << "NaN";
                  return os;
               }

               if (latex_as_equation())
               {
                  if (!inside_equation)
                     os.put(os.widen('$'));
                  os << mantissa;
                  if (exponent.size())
                  {
                     if (this->multiply() == multiply_x)
                        os << "\\mathrm{x}";
                     else if (this->multiply() == multiply_dot)
                        os << "\\cdot";
                     else
                        os << "\\times";
                     os << " 10^{" << exponent << "}";
                  }
                  if (!inside_equation)
                     os.put(os.widen('$'));
               }
               else
               {
                  os << mantissa;
                  if (exponent.size())
                  {
                     if (this->multiply() == multiply_x)
                        os << "x";
                     else if (this->multiply() == multiply_dot)
                        os << "$\\cdot$";
                     else
                        os << "$\\times$";
                     os << "10\\textsuperscript{" << exponent << "}";
                  }
               }
               return os;
            }
            template <class Complex>
            std::basic_ostream<charT, Traits>& format_complex(std::basic_ostream<charT, Traits>& os, const Complex& f)
            {
               if (!inside_equation && latex_as_equation())
                  os.put(os.widen('$'));
               bool saved_inside_equation = inside_equation;
               inside_equation = latex_as_equation();
               if ((std::isnan(f.real())) || (std::isnan)(f.imag()))
               {
                  if (inside_equation)
                     os << "\\mathrm{NaN}";
                  else
                     os << "NaN";
                  return os;
               }
               else if ((std::isinf)(f.real()) || (std::isinf)(f.imag()))
               {
                  if (!inside_equation)
                     os.put(os.widen('$'));
                  os << "\\tilde{\\infty}";
                  if (!inside_equation)
                     os.put(os.widen('$'));
                  return os;
               }
               else
               {
                  bool need_i = true;

                  if (!this->show_zero_components() && (f.imag() == 0))
                  {
                     format_float(os, f.real());
                     need_i = false;
                  }
                  else if (!this->show_zero_components() && (f.real() == 0))
                  {
                     format_float(os, f.imag());
                  }
                  else
                  {

                     format_float(os, f.real());
                     typename Complex::value_type i(f.imag());
                     bool isneg = i < 0;
                     if (isneg)
                     {
                        i = -i;
                        os << " - ";
                     }
                     else
                        os << " + ";
                     format_float(os, i);
                  }
                  if (need_i)
                  {
                     if (this->imaginary_style() == doublestruck_i)
                        BOOST_THROW_EXCEPTION(std::runtime_error("Doublestruck imaginary i not implemented for LaTex output."));
                     else if ((this->imaginary_style() == upright_i) && latex_as_equation())
                        os << "\\mathrm{i}";
                     else if ((this->imaginary_style() == slanted_i) && !latex_as_equation())
                        os << "\\textit{i}";
                     else
                        os.put(os.widen('i'));
                  }
               }
               inside_equation = saved_inside_equation;
               if (!inside_equation && latex_as_equation())
                  os.put(os.widen('$'));

               return os;
            }
         };

         template <class charT, class Traits>
         class basic_numeric_formatter<html_format, charT, Traits> : public basic_numeric_formatter_base<charT, Traits>
         {
         protected:
            template <class Integer>
            std::basic_ostream<charT, Traits>& format_integer(std::basic_ostream<charT, Traits>& os, const Integer& i)
            {
               if (this->styling() >= minimal_styling)
                  os << "<span class=\"number\">";
               if (this->styling() > minimal_styling)
                  os << "<span class=\"integer\">";
               os << i;
               if (this->styling() > minimal_styling)
                  os << "</span>";
               if (this->styling() >= minimal_styling)
                  os << "</span>";
               return os;
            }
            template <class Float>
            std::basic_ostream<charT, Traits>& format_float(std::basic_ostream<charT, Traits>& os, const Float& f)
            {
               if (this->styling() >= minimal_styling)
                  os << "<span class=\"number\">";
               if (this->styling() > minimal_styling)
                  os << "<span class=\"float\">";

               if ((std::isinf)(f))
               {
                  if (f < 0)
                     os.put(os.widen('-'));
                  os << "&#8734;";
               }
               else if ((std::isnan(f)))
               {
                  os << "NaN";
               }
               else
               {
                  std::basic_string<charT, Traits> mantissa, exponent;
                  basic_numeric_formatter_base<charT, Traits>::decompose_float(os, f, mantissa, exponent);

                  if (exponent.size())
                  {
                     os << mantissa;
                     if (this->multiply() == multiply_x)
                        os.put(os.widen('x'));
                     else if (this->multiply() == multiply_dot)
                        os << "&#x22C5;";
                     else
                        os << "&#x00D7;";
                     os << "10<sup>" << exponent << "</sup>";
                  }
                  else
                     os << mantissa;
               }

               if (this->styling() > minimal_styling)
                  os << "</span>";
               if (this->styling() >= minimal_styling)
                  os << "</span>";

               return os;
            }
            template <class Complex>
            std::basic_ostream<charT, Traits>& format_complex(std::basic_ostream<charT, Traits>& os, const Complex& f)
            {
               if (this->styling() >= minimal_styling)
                  os << "<span class=\"number\">";
               if (this->styling() > minimal_styling)
                  os << "<span class=\"complex\">";

               if ((std::isnan(f.real())) || (std::isnan)(f.imag()))
               {
                  os << "NaN";
               }
               else if ((std::isinf)(f.real()) || (std::isinf)(f.imag()))
               {
                  os << "&#8734;&#x0303;";
               }
               else
               {
                  styling_level_t saved_style = this->styling();
                  this->styling(no_styling);
                  
                  bool need_i = true;

                  if(!this->show_zero_components() && (f.imag() == 0))
                  { 
                     format_float(os, f.real());
                     need_i = false;
                  }
                  else if (!this->show_zero_components() && (f.real() == 0))
                  {
                     format_float(os, f.imag());
                  }
                  else
                  {
                     format_float(os, f.real());
                     typename Complex::value_type i(f.imag());
                     bool isneg = i < 0;
                     if (isneg)
                     {
                        i = -i;
                        os << " - ";
                     }
                     else
                        os << " + ";
                     format_float(os, i);
                  }
                  if (need_i)
                  {
                     if (saved_style >= minimal_styling)
                        os << "<span class=\"imaginary_i\">";
                     if (this->imaginary_style() == doublestruck_i)
                        os << "&#x2148;";
                     else if (this->imaginary_style() == slanted_i)
                        os << "<i>i</i>";
                     else
                        os.put(os.widen('i'));
                     if (saved_style >= minimal_styling)
                        os << "</span>";
                  }
                  this->styling(saved_style);
               }

               if (this->styling() > minimal_styling)
                  os << "</span>";
               if (this->styling() >= minimal_styling)
                  os << "</span>";

               return os;
            }
         };


         template <output_format_t Format = text_format, class charT = char, class Traits = std::char_traits<charT> >
         class basic_numeric_printer : private basic_numeric_formatter<Format, charT, Traits>
         {
            std::basic_ostream<charT, Traits>* p_stream;
            typedef basic_numeric_formatter<Format, charT, Traits> base_type;

         public:
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os) : p_stream(&os) {}
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os, const basic_numeric_printer& printer) 
               : p_stream(&os) 
            {
               styling(printer.styling());
               imaginary_style(printer.imaginary_style());
               latex_as_equation(printer.latex_as_equation());
               parenthesis(printer.parenthesis());
               multiply(printer.multiply());
               show_zero_components(printer.show_zero_components());
            }

            using base_type::format_integer;
            using base_type::format_float;
            using base_type::format_complex;
            using base_type::styling;
            using base_type::imaginary_style;
            using base_type::latex_as_equation;
            using base_type::parenthesis;
            using base_type::multiply;
            using base_type::show_zero_components;

            std::basic_ostream<charT, Traits>& stream()
            {
               return *p_stream;
            }
         };

         typedef basic_numeric_printer<> text_printer;
         typedef basic_numeric_printer<docbook_format> docbook_printer;
         typedef basic_numeric_printer<latex_format> latex_printer;
         typedef basic_numeric_printer<html_format> html_printer;

         typedef basic_numeric_printer<text_format, wchar_t> wtext_printer;
         typedef basic_numeric_printer<docbook_format, wchar_t> wdocbook_printer;
         typedef basic_numeric_printer<latex_format, wchar_t> wlatex_printer;
         typedef basic_numeric_printer<html_format, wchar_t> whtml_printer;

         template <output_format_t Format, class charT, class Traits, class Integer>
         typename boost::enable_if_c<is_integer_like<Integer>::value, basic_numeric_printer<Format, charT, Traits>&>::type operator << (basic_numeric_printer<Format, charT, Traits>& os, const Integer& i)
         {
            if (os.stream().width())
            {
               std::basic_stringstream<charT, Traits> ss;
               ss.copyfmt(os.stream());
               ss.imbue(os.stream().getloc());
               ss.width(0);
               basic_numeric_printer<Format, charT, Traits> fmt(ss, os);
               fmt.format_integer(fmt.stream(), i);
               os.stream() << ss.str();
            }
            else
               os.format_integer(os.stream(), i);
            return os;
         }

         template <output_format_t Format, class charT, class Traits, class Float>
         typename boost::enable_if_c<is_float_like<Float>::value, basic_numeric_printer<Format, charT, Traits>&>::type operator << (basic_numeric_printer<Format, charT, Traits>& os, const Float& f)
         {
            if (os.stream().width())
            {
               std::basic_stringstream<charT, Traits> ss;
               ss.copyfmt(os.stream());
               ss.imbue(os.stream().getloc());
               ss.width(0);
               basic_numeric_printer<Format, charT, Traits> fmt(ss, os);
               fmt.format_float(fmt.stream(), f);
               os.stream() << ss.str();
            }
            else
               os.format_float(os.stream(), f);
            return os;
         }

         template <output_format_t Format, class charT, class Traits, class Complex>
         typename boost::enable_if_c<is_complex_like<Complex>::value, basic_numeric_printer<Format, charT, Traits>&>::type operator << (basic_numeric_printer<Format, charT, Traits>& os, const Complex& f)
         {
            if (os.stream().width())
            {
               std::basic_stringstream<charT, Traits> ss;
               ss.copyfmt(os.stream());
               ss.imbue(os.stream().getloc());
               ss.width(0);
               basic_numeric_printer<Format, charT, Traits> fmt(ss, os);
               fmt.format_complex(fmt.stream(), f);
               os.stream() << ss.str();
            }
            else
               os.format_complex(os.stream(), f);
            return os;
         }

         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, styling_level_t s)
         {
            os.styling(s);
            return os;
         }

         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, multiplyer_t s)
         {
            os.multiply(s);
            return os;
         }

         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, const latex_as_equation_t* const&)
         {
            os.latex_as_equation(true);
            return os;
         }

         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, const latex_as_text_t* const&)
         {
            os.latex_as_equation(false);
            return os;
         }

         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, const imaginary_i_t i)
         {
            os.imaginary_style(i);
            return os;
         }

         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, const zero_component_t z)
         {
            os.show_zero_components(z == ::boost::math::tools::show_zero_components ? true : false);
            return os;
         }

         //
         // Output formatters:
         //
         template <output_format_t Format, class charT, class Traits, class T>
         typename boost::enable_if_c<is_unhandled_type<T>::value, basic_numeric_printer<Format, charT, Traits>&>::type operator<< (basic_numeric_printer<Format, charT, Traits>& os, const T& unhandled)
         {
            os.stream() << unhandled;
            return os;
         }
         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator<< (basic_numeric_printer<Format, charT, Traits>& os, std::basic_ostream<charT, Traits>& (*pf)(std::basic_ostream<charT, Traits>&))
         {
            pf(os.stream());
            return os;
         }
         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator<< (basic_numeric_printer<Format, charT, Traits>& os, std::basic_ios<charT, Traits>& (*pf)(std::basic_ios<charT, Traits>&))
         {
            pf(os.stream());
            return os;
         }
         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator<< (basic_numeric_printer<Format, charT, Traits>& os, std::ios_base& (*pf)(std::ios_base&))
         {
            pf(os.stream());
            return os;
         }
      }
   }
}

#endif // BOOST_MATH_TOOLS_FORMATTING_HPP

