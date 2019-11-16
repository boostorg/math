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
#include <boost/regex/pending/unicode_iterator.hpp>

#ifdef _MSC_VER
#pragma once
#endif

namespace boost {

   template <class T>
   struct rational;

   namespace math {
      namespace tools {

         template <class T>
         T numerator(const boost::rational<T>& r) { return r.numerator(); }
         template <class T>
         T denominator(const boost::rational<T>& r) { return r.denominator(); }

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
         using polynomial_t = decltype(std::declval<const T&>()(std::declval<const T&>()[0]));
         template <class T>
         using polynomial2_t = decltype(std::declval<const T&>().degree());
         template <class T>
         using rational_t = decltype(numerator(std::declval<const T&>()) + denominator(std::declval<const T&>()));

         template <class T>
         struct is_complex_like
         {
            static const bool value = boost::is_detected_v<real_t, T> && boost::is_detected_v<imag_t, T>;
         };

         template <class T>
         struct is_polynomial_like
         {
            static const bool value = boost::is_detected_v<polynomial_t, T> && boost::is_detected_v<polynomial2_t, T>;
         };

         template <class T>
         struct is_rational_like
         {
            static const bool value = boost::is_detected_v<rational_t, T>;
         };

         template <class T>
         struct is_unhandled_type
         {
            static const bool value = !is_complex_like<T>::value && !is_float_like<T>::value && !is_integer_like<T>::value && !is_polynomial_like<T>::value && !is_rational_like<T>::value;
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

         enum latex_t
         {
            latex_as_text = 1,
            latex_as_equation = 0
         };

         enum unicode_text_t
         {
            unicode_text_output = 0,
            ascii_text_output = 1
         };

         template<class charT, class Traits>
         typename boost::enable_if_c<sizeof(charT) == 1>::type write_unicode_char(std::basic_ostream <charT, Traits>& os, std::uint32_t code_point)
         {
            std::uint32_t u[2] = { code_point, 0 };
            boost::u32_to_u8_iterator<const std::uint32_t*, charT> i(u), j(&u[0] + 1);
            while (i != j)
            {
               os.put(*i++);
            }
         }
         template<class charT, class Traits>
         typename boost::enable_if_c<sizeof(charT) == 2>::type write_unicode_char(std::basic_ostream <charT, Traits>& os, std::uint32_t code_point)
         {
            std::uint32_t u[2] = { code_point, 0 };
            boost::u32_to_u16_iterator<const std::uint32_t*, charT> i(u), j(&u[0] + 1);
            while (i != j)
            {
               os.put(*i++);
            }
         }
         template<class charT, class Traits>
         typename boost::enable_if_c<sizeof(charT) == 4>::type write_unicode_char(std::basic_ostream <charT, Traits>& os, std::uint32_t code_point)
         {
            os.put(code_point);
         }

         template<class charT, class Traits>
         void write_unicode_superscript(std::basic_ostream <charT, Traits>& os, std::basic_string<charT, Traits> digits)
         {
            for (std::size_t i = 0; i < digits.size(); ++i)
            {
               switch (digits[i])
               {
               case '0':
                  write_unicode_char(os, 0x2070);
                  break;
               case '1':
                  write_unicode_char(os, 0xB9);
                  break;
               case '2':
                  write_unicode_char(os, 0xB2);
                  break;
               case '3':
                  write_unicode_char(os, 0xB3);
                  break;
               case '4':
                  write_unicode_char(os, 0x2074);
                  break;
               case '5':
                  write_unicode_char(os, 0x2075);
                  break;
               case '6':
                  write_unicode_char(os, 0x2076);
                  break;
               case '7':
                  write_unicode_char(os, 0x2076);
                  break;
               case '8':
                  write_unicode_char(os, 0x2078);
                  break;
               case '9':
                  write_unicode_char(os, 0x2079);
                  break;
               case '+':
                  write_unicode_char(os, 0x207A);
                  break;
               case '-':
                  write_unicode_char(os, 0x207B);
                  break;
               }
            }
         }

         template<class charT, class Traits>
         void write_unicode_subscript(std::basic_ostream <charT, Traits>& os, std::basic_string<charT, Traits> digits)
         {
            for (std::size_t i = 0; i < digits.size(); ++i)
            {
               if ((digits[i] >= '0') && (digits[i] <= '9'))
                  write_unicode_char(os, 0x2080 + (digits[i] - '0'));
               else if (digits[i] == '-')
                  write_unicode_char(os, 0x208B);
               else if (digits[i] == '+')
                  write_unicode_char(os, 0x208A);
               else
                  os.put(digits[i]);
            }
         }

         template <class charT>
         typename boost::enable_if_c<sizeof(charT) == 1, std::size_t>::type unicode_character_len(const charT* p)
         {
            std::size_t l = 0;
            for (unsigned i = 0; p && p[i]; ++i) 
               ++l;
            boost::u8_to_u32_iterator<const charT*> i(p), j(p + l);
            return std::distance(i, j);
         }
         template <class charT>
         typename boost::enable_if_c<sizeof(charT) == 2, std::size_t>::type unicode_character_len(const charT* p)
         {
            std::size_t l = 0;
            for (unsigned i = 0; p && p[i]; ++i) 
               ++l;
            boost::u16_to_u32_iterator<const charT*> i(p), j(p + l);
            return std::distance(i, j);
         }
         template <class charT>
         typename boost::enable_if_c<sizeof(charT) == 4, std::size_t>::type unicode_character_len(const charT* p)
         {
            std::size_t l = 0;
            for (unsigned i = 0; p && p[i]; ++i) 
               ++l;
            return l;
         }

         template <class T>
         typename boost::enable_if_c<is_complex_like<T>::value, bool>::type iszero(const T& val)
         {
            return (val.real() == 0) && (val.imag() == 0);
         }
         template <class T>
         typename boost::disable_if_c<is_complex_like<T>::value, bool>::type iszero(const T& val)
         {
            return val == 0;
         }
         template <class T>
         typename boost::enable_if_c<is_complex_like<T>::value, bool>::type isneg(const T& val)
         {
            return (val.real() < 0);
         }
         template <class T>
         typename boost::disable_if_c<is_complex_like<T>::value, bool>::type isneg(const T& val)
         {
            return val < 0;
         }

         template <output_format_t Format, class charT, class Traits>
         class basic_numeric_formatter;

         template <class charT, class Traits>
         class basic_numeric_formatter_base
         {
            styling_level_t styling_level;
            imaginary_i_t i_style;
            multiplyer_t  multiply_style;
            std::size_t m_parenthesis;
            bool m_show_zero_components;
            bool m_use_unicode;
            std::basic_ostream<charT, Traits>* p_stream;
         protected:
            template <class Float>
            static void decompose_float(std::basic_ostream<charT, Traits>& os, const Float& f, std::basic_string<charT, Traits>& mantissa, std::basic_string<charT, Traits>& exponent)
            {
               static const charT exponent_string[3] = { 'e', 'E', 0 };
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
            basic_numeric_formatter_base(std::basic_ostream<charT, Traits>& os) 
               : styling_level(full_styling), i_style(upright_i), multiply_style(multiply_times), 
                 m_parenthesis(0), m_show_zero_components(false), m_use_unicode(true), p_stream(&os) {}
            basic_numeric_formatter_base(std::basic_ostream<charT, Traits>& os, const basic_numeric_formatter_base& o) 
               : styling_level(o.styling_level), i_style(o.i_style), multiply_style(o.multiply_style),
                 m_parenthesis(o.m_parenthesis), m_show_zero_components(o.m_show_zero_components), m_use_unicode(o.m_use_unicode), p_stream(&os) {}

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
            void parenthesis(std::size_t b)
            {
               m_parenthesis = b;
            }
            std::size_t parenthesis()const
            {
               return m_parenthesis;
            }
            void show_zero_components(bool b)
            {
               m_show_zero_components = b;
            }
            bool show_zero_components()const
            {
               return m_show_zero_components;
            }
            void use_unicode(bool b)
            {
               m_use_unicode = b;
            }
            bool use_unicode()const
            {
               return m_use_unicode;
            }

            std::basic_ostream<charT, Traits>& stream()
            {
               return *p_stream;
            }

            struct scoped_parenthesis
            {
               basic_numeric_formatter_base<charT, Traits>* m_formatter;
               scoped_parenthesis(basic_numeric_formatter_base<charT, Traits>* formatter) : m_formatter(formatter)
               {
                  m_formatter->parenthesis(1 + m_formatter->parenthesis());
               }
               ~scoped_parenthesis()
               {
                  m_formatter->parenthesis(m_formatter->parenthesis() - 1);
               }
            };

            template <class Integer>
            static std::basic_ostream<charT, Traits>& format_integer(std::basic_ostream<charT, Traits>& os, const Integer& i)
            {
               return os << i;
            }
            template <class Float>
            std::basic_ostream<charT, Traits>& format_float(std::basic_ostream<charT, Traits>& os, const Float& f)
            {
               if ((std::isinf)(f))
               {
                  if (f < 0)
                     os << "-";
                  if (use_unicode())
                  {
                     write_unicode_char(os, 0x221E);
                  }
                  else
                     os << "INFINITY";
                  return os;
               }
               if ((std::isnan)(f))
               {
                  return os << "NaN";
               }
               std::basic_string<charT, Traits> mantissa, exponent;
               decompose_float(os, f, mantissa, exponent);
               os << mantissa;
               if (exponent.size())
               {
                  if (use_unicode())
                  {
                     if (multiply() == multiply_times)
                        write_unicode_char(os, 0xD7);
                     else if (multiply() == multiply_dot)
                        write_unicode_char(os, 0x22C5);
                     else
                        os << "x";
                     os << "10";
                     write_unicode_superscript(os, exponent);
                  }
                  else
                     os << "x10^" << exponent;
               }
               return os;
            }
            template <class Complex>
            std::basic_ostream<charT, Traits>& format_complex(std::basic_ostream<charT, Traits>& os, const Complex& f)
            {
               if ((std::isnan)(f.real()) || (std::isnan)(f.imag()))
               {
                  return os << "NaN";
               }
               if ((std::isinf)(f.real()) || (std::isinf)(f.imag()))
               {
                  if (use_unicode())
                  {
                     write_unicode_char(os, 0x221E);
                     write_unicode_char(os, 0x0303);
                  }
                  else
                     os << "COMPLEX INFINITY";
                  return os;
               }
               if ((!this->show_zero_components()) && (f.imag() == 0))
               {
                  format_float(os, f.real());
                  return os;
               }
               else if ((!this->show_zero_components()) && (f.real() == 0))
               {
                  format_float(os, f.imag());
               }
               else
               {
                  if (parenthesis())
                     os << "(";
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
               if ((use_unicode()) && (this->imaginary_style() == doublestruck_i))
                  write_unicode_char(os, 0x2148);
               else if ((use_unicode()) && (this->imaginary_style() == slanted_i))
                  write_unicode_char(os, 0x1D456);
               else
                  os.put(os.widen('i'));
               if (parenthesis() && !((!this->show_zero_components()) && (f.real() == 0)))
                  os << ")";
               return os;
            }
            template <class Polynomial>
            std::basic_ostream<charT, Traits>& format_polynomial(std::basic_ostream<charT, Traits>& os, const Polynomial& f)
            {
               scoped_parenthesis scoped(this);
               bool have_first = false;
               for (unsigned i = 0; i <= f.degree(); ++i)
               {
                  auto coef = f[i];
                  if (show_zero_components() || (!iszero(coef)))
                  {
                     if (have_first)
                     {
                        if (!isneg(coef))
                           os << " + ";
                        else
                           os << " - ";
                     }
                     print(*this, os, !isneg(coef) ? coef : -coef);
                     have_first = true;
                     if (i)
                     {
                        os << "x";
                        if (i > 1)
                        {
                           if (use_unicode())
                           {
                              std::basic_stringstream<charT, Traits> ss;
                              ss << i;
                              std::basic_string<charT, Traits> s = ss.str();
                              write_unicode_superscript(os, s);
                           }
                           else
                              os << "^" << i;
                        }
                     }
                  }
               }
               return os;
            }

            template <class Rational>
            std::basic_ostream<charT, Traits>& format_rational(std::basic_ostream<charT, Traits>& os, const Rational& rat)
            {
               auto num = numerator(rat);
               auto denom = denominator(rat);


               if (num == 0)
               {
                  return os << "0";
               }
               else if (denom == 1)
               {
                  format_integer(os, num);
                  return os;
               }
               else
               {
                  if (parenthesis())
                     os << "(";
                  auto s1 = part_as_string(*this, num);
                  auto s2 = part_as_string(*this, denom);

                  write_unicode_superscript(os, s1);
                  os << "/";
                  write_unicode_subscript(os, s2);

                  if (parenthesis())
                     os << ")";
               }

               return os;
            }

            bool latex_as_equation()const { return false; }
            void latex_as_equation(bool) {}

            template <class Printer, class Value>
            static std::basic_string<charT, Traits> part_as_string(Printer& printer, const Value& value)
            {
               std::basic_ostringstream<charT, Traits> ss;
               ss.copyfmt(printer.stream());
               ss.imbue(printer.stream().getloc());
               ss.width(0);
               Printer fmt(ss, printer);
               fmt.print(fmt, fmt.stream(), value);
               return ss.str();
            }

            template <class Printer, class Value>
            static typename boost::enable_if_c<is_integer_like<Value>::value>::type print(Printer& printer, std::basic_ostream<charT, Traits>& os, const Value& value)
            {
               printer.format_integer(os, value);
            }
            template <class Printer, class Value>
            static typename boost::enable_if_c<is_float_like<Value>::value>::type print(Printer& printer, std::basic_ostream<charT, Traits>& os, const Value& value)
            {
               printer.format_float(os, value);
            }
            template <class Printer, class Value>
            static typename boost::enable_if_c<is_complex_like<Value>::value>::type print(Printer& printer, std::basic_ostream<charT, Traits>& os, const Value& value)
            {
               printer.format_complex(os, value);
            }
            template <class Printer, class Value>
            static typename boost::enable_if_c<is_polynomial_like<Value>::value>::type print(Printer& printer, std::basic_ostream<charT, Traits>& os, const Value& value)
            {
               printer.format_polynomial(os, value);
            }
            template <class Printer, class Value>
            static typename boost::enable_if_c<is_rational_like<Value>::value>::type print(Printer& printer, std::basic_ostream<charT, Traits>& os, const Value& value)
            {
               printer.format_rational(os, value);
            }
         };

         template <class charT, class Traits>
         class basic_numeric_formatter<text_format, charT, Traits> : public basic_numeric_formatter_base<charT, Traits>
         {
         public:
            basic_numeric_formatter(std::basic_ostream<charT, Traits>& os) : basic_numeric_formatter_base<charT, Traits>(os) {}
         };

         template <class charT, class Traits>
         class basic_numeric_formatter<docbook_format, charT, Traits> : public basic_numeric_formatter_base<charT, Traits>
         {
            friend class basic_numeric_formatter_base<charT, Traits>;
         protected:
            basic_numeric_formatter(std::basic_ostream<charT, Traits>& os) : basic_numeric_formatter_base<charT, Traits>(os) {}

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
                  bool need_paren = false;

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
                     if (this->parenthesis())
                     {
                        os << "(";
                        need_paren = true;
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
                  if (need_paren)
                     os << ")";
                  this->styling(saved_style);
               }

               if (this->styling() > minimal_styling)
                  os << "</phrase>";
               if (this->styling() >= minimal_styling)
                  os << "</phrase>";

               return os;
            }
            template <class Polynomial>
            std::basic_ostream<charT, Traits>& format_polynomial(std::basic_ostream<charT, Traits>& os, const Polynomial& f)
            {
               typename basic_numeric_formatter_base<charT, Traits>::scoped_parenthesis scoped(this);

               if (this->styling() >= minimal_styling)
                  os << "<phrase role=\"number\">";
               if (this->styling() > minimal_styling)
                  os << "<phrase role=\"polynomial\">";

               bool have_first = false;
               for (unsigned i = 0; i <= f.degree(); ++i)
               {
                  auto coef = f[i];
                  if (this->show_zero_components() || !iszero(coef))
                  {
                     if (have_first)
                     {
                        if (!isneg(coef))
                           os << " + ";
                        else
                           os << " - ";
                     }
                     basic_numeric_formatter_base<charT, Traits>::print(*this, os, !isneg(coef) ? coef : -coef);
                     have_first = true;
                     if (i)
                     {
                        os << "<emphasis>x</emphasis>";
                        if (i > 1)
                        {
                           os << "<superscript>" << i << "</superscript>";
                        }
                     }
                  }
               }

               if (this->styling() > minimal_styling)
                  os << "</phrase>";
               if (this->styling() >= minimal_styling)
                  os << "</phrase>";

               return os;
            }
            template <class Rational>
            std::basic_ostream<charT, Traits>& format_rational(std::basic_ostream<charT, Traits>& os, const Rational& rat)
            {
               auto num = numerator(rat);
               auto denom = denominator(rat);

               if ((denom == 1) || (num == 0))
               {
                  format_integer(os, num);
               }
               else
               {
                  if (this->styling() >= minimal_styling)
                     os << "<phrase role=\"number\">";
                  if (this->styling() > minimal_styling)
                     os << "<phrase role=\"rational\">";
                  
                  if (this->parenthesis())
                     os << "(";
                  
                  os << "<superscript>" << num << "</superscript>";
                  os << "&#x2044;<subscript>" << denom << "</subscript>";

                  if (this->parenthesis())
                     os << ")";

                  if (this->styling() > minimal_styling)
                     os << "</phrase>";
                  if (this->styling() >= minimal_styling)
                     os << "</phrase>";
               }
               return os;
            }
         };

         template <class charT, class Traits>
         class basic_numeric_formatter<latex_format, charT, Traits> : public basic_numeric_formatter_base<charT, Traits>
         {
            friend class basic_numeric_formatter_base<charT, Traits>;

            bool is_latex_as_equation;
            bool inside_equation;

         public:
            basic_numeric_formatter(std::basic_ostream<charT, Traits>& os) : basic_numeric_formatter_base<charT, Traits>(os), is_latex_as_equation(true), inside_equation(false) {}
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
                  bool need_paren = false;

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
                     if (this->parenthesis())
                     {
                        os << "(";
                        need_paren = true;
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
                  if (need_paren)
                     os << ")";
               }
               inside_equation = saved_inside_equation;
               if (!inside_equation && latex_as_equation())
                  os.put(os.widen('$'));

               return os;
            }
            template <class Polynomial>
            std::basic_ostream<charT, Traits>& format_polynomial(std::basic_ostream<charT, Traits>& os, const Polynomial& f)
            {
               typename basic_numeric_formatter_base<charT, Traits>::scoped_parenthesis scoped(this);

               if (!inside_equation && latex_as_equation())
                  os.put(os.widen('$'));
               bool saved_inside_equation = inside_equation;
               inside_equation = latex_as_equation();
               bool have_first = false;

               for (unsigned i = 0; i <= f.degree(); ++i)
               {
                  auto coef = f[i];
                  if (this->show_zero_components() || !iszero(coef))
                  {
                     if (have_first)
                     {
                        if (!isneg(coef))
                           os << " + ";
                        else
                           os << " - ";
                     }
                     basic_numeric_formatter_base<charT, Traits>::print(*this, os, !isneg(coef) ? coef : -coef);
                     have_first = true;
                     if (i)
                     {
                        if (inside_equation)
                           os << "x";
                        else
                           os << "\\textit{x}";
                        if (i > 1)
                        {
                           if (inside_equation)
                              os << "^{" << i << "}";
                           else
                              os << "\\textsuperscript{" << i << "}";
                        }
                     }
                  }
               }

               inside_equation = saved_inside_equation;
               if (!inside_equation && latex_as_equation())
                  os.put(os.widen('$'));

               return os;
            }
            template <class Rational>
            std::basic_ostream<charT, Traits>& format_rational(std::basic_ostream<charT, Traits>& os, const Rational& rat)
            {
               auto num = numerator(rat);
               auto denom = denominator(rat);


               if ((num == 0) || (denom == 1))
               {
                  format_integer(os, num);
               }
               else
               {
                  if (!inside_equation && latex_as_equation())
                     os.put(os.widen('$'));
                  if (this->parenthesis())
                     os << "(";
                  
                  if (latex_as_equation())
                  {
                     if(num < 0)
                        os << "-\\frac{" << -num << "}{" << denom << "}";
                     else
                        os << "\\frac{" << num << "}{" << denom << "}";
                  }
                  else
                     os << "\\textsuperscript{ " << num << "}/\\textsubscript{" << denom << "}";

                  if (this->parenthesis())
                     os << ")";
                  if (!inside_equation && latex_as_equation())
                     os.put(os.widen('$'));
               }

               return os;
            }
         };

         template <class charT, class Traits>
         class basic_numeric_formatter<html_format, charT, Traits> : public basic_numeric_formatter_base<charT, Traits>
         {
            friend class basic_numeric_formatter_base<charT, Traits>;
         protected:
            basic_numeric_formatter(std::basic_ostream<charT, Traits>& os) : basic_numeric_formatter_base<charT, Traits>(os) {}

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
                  bool need_paren = false;

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
                     if (this->parenthesis())
                     {
                        os << "(";
                        need_paren = true;
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
                  if (need_paren)
                     os << ")";
                  this->styling(saved_style);
               }

               if (this->styling() > minimal_styling)
                  os << "</span>";
               if (this->styling() >= minimal_styling)
                  os << "</span>";

               return os;
            }
            template <class Polynomial>
            std::basic_ostream<charT, Traits>& format_polynomial(std::basic_ostream<charT, Traits>& os, const Polynomial& f)
            {
               typename basic_numeric_formatter_base<charT, Traits>::scoped_parenthesis scoped(this);

               if (this->styling() >= minimal_styling)
                  os << "<span class=\"number\">";
               if (this->styling() > minimal_styling)
                  os << "<span class=\"polynomial\">";

               styling_level_t saved_style = this->styling();
               this->styling(no_styling);

               bool have_first = false;
               for (unsigned i = 0; i <= f.degree(); ++i)
               {
                  auto coef = f[i];
                  if (this->show_zero_components() || !iszero(coef))
                  {
                     if (have_first)
                     {
                        if (!isneg(coef))
                           os << " + ";
                        else
                           os << " - ";
                     }
                     basic_numeric_formatter_base<charT, Traits>::print(*this, os, !isneg(coef) ? coef : -coef);
                     have_first = true;
                     if (i)
                     {
                        os << "<i>x</i>";
                        if (i > 1)
                        {
                           os << "<sup>" << i << "</sup>";
                        }
                     }
                  }
               }

               this->styling(saved_style);

               if (this->styling() > minimal_styling)
                  os << "</span>";
               if (this->styling() >= minimal_styling)
                  os << "</span>";

               return os;
            }
            template <class Rational>
            std::basic_ostream<charT, Traits>& format_rational(std::basic_ostream<charT, Traits>& os, const Rational& rat)
            {
               auto num = numerator(rat);
               auto denom = denominator(rat);

               if ((num == 0) || (denom == 1))
               {
                  format_integer(os, num);
               }
               else
               {
                  if (this->styling() >= minimal_styling)
                     os << "<span class=\"number\">";
                  if (this->styling() > minimal_styling)
                     os << "<span class=\"rational\">";

                  if (this->parenthesis())
                     os << "(";
                  os << "<sup>" << num << "</sup>";
                  os << "&#x2044;<sub>" << denom << "</sub>";

                  if (this->parenthesis())
                     os << ")";

                  if (this->styling() > minimal_styling)
                     os << "</span>";
                  if (this->styling() >= minimal_styling)
                     os << "</span>";
               }

               return os;
            }
         };


         template <output_format_t Format = text_format, class charT = char, class Traits = std::char_traits<charT> >
         class basic_numeric_printer : private basic_numeric_formatter<Format, charT, Traits>
         {
            typedef basic_numeric_formatter<Format, charT, Traits> base_type;

         public:
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os) : basic_numeric_formatter<Format, charT, Traits>(os) {}
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os, const basic_numeric_printer& printer) 
               : basic_numeric_formatter<Format, charT, Traits>(os)
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
            using base_type::format_polynomial;
            using base_type::format_rational;
            using base_type::styling;
            using base_type::imaginary_style;
            using base_type::latex_as_equation;
            using base_type::parenthesis;
            using base_type::multiply;
            using base_type::show_zero_components;
            using base_type::use_unicode;
            using base_type::print;
            using base_type::stream;
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
         typename boost::enable_if_c<!is_unhandled_type<Integer>::value, basic_numeric_printer<Format, charT, Traits>&>::type operator << (basic_numeric_printer<Format, charT, Traits>& os, const Integer& i)
         {
            std::size_t w = (std::size_t)os.stream().width();
            if (w)
            {
               std::basic_string<charT, Traits> s = basic_numeric_formatter_base<charT, Traits>::part_as_string(os, i);
               std::size_t len = unicode_character_len(s.c_str());
               if (len < w)
               {
                  auto pos = os.stream().flags() & std::ios_base::adjustfield;
                  s.insert(pos == std::ios_base::left ? s.end() : s.begin(), w - len, os.stream().fill());
               }
               os.stream() << s;
            }
            else
               os.print(os, os.stream(), i);
            return os;
         }
#if 0
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

         template <output_format_t Format, class charT, class Traits, class Polynomial>
         typename boost::enable_if_c<is_polynomial_like<Polynomial>::value, basic_numeric_printer<Format, charT, Traits>&>::type operator << (basic_numeric_printer<Format, charT, Traits>& os, const Polynomial& f)
         {
            if (os.stream().width())
            {
               std::basic_stringstream<charT, Traits> ss;
               ss.copyfmt(os.stream());
               ss.imbue(os.stream().getloc());
               ss.width(0);
               basic_numeric_printer<Format, charT, Traits> fmt(ss, os);
               fmt.format_polynomial(fmt.stream(), f);
               os.stream() << ss.str();
            }
            else
               os.format_polynomial(os.stream(), f);
            return os;
         }
#endif
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
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, latex_t l)
         {
            os.latex_as_equation(l == latex_as_equation);
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

         template <output_format_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, const unicode_text_t z)
         {
            os.use_unicode(z == ::boost::math::tools::unicode_text_output ? true : false);
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

