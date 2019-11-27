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

         namespace detail {

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
            using interval_t = decltype(singleton(std::declval<const T&>()) + intersect(std::declval<const T&>(), std::declval<const T&>()) + median(std::declval<const T&>()));

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
            struct is_interval_like
            {
               static const bool value = boost::is_detected_v<interval_t, T>;
            };

            template <class T>
            struct is_float_like<T, true>
            {
               // The point "floats":
               static const bool value = (std::numeric_limits<T>::max_exponent != std::numeric_limits<T>::min_exponent) && !is_interval_like<T>::value;
            };

            template <class T>
            struct is_unhandled_type
            {
               static const bool value = !is_complex_like<T>::value && !is_float_like<T>::value && !is_integer_like<T>::value && !is_polynomial_like<T>::value && !is_rational_like<T>::value;
            };

            template<class charT, class Traits>
            typename boost::enable_if_c<sizeof(charT) == 1>::type print_unicode_char(std::basic_ostream <charT, Traits>& os, std::uint32_t code_point)
            {
               std::uint32_t u[2] = { code_point, 0 };
               boost::u32_to_u8_iterator<const std::uint32_t*, charT> i(u), j(&u[0] + 1);
               while (i != j)
               {
                  os.put(*i++);
               }
            }
            template<class charT, class Traits>
            typename boost::enable_if_c<sizeof(charT) == 2>::type print_unicode_char(std::basic_ostream <charT, Traits>& os, std::uint32_t code_point)
            {
               std::uint32_t u[2] = { code_point, 0 };
               boost::u32_to_u16_iterator<const std::uint32_t*, charT> i(u), j(&u[0] + 1);
               while (i != j)
               {
                  os.put(*i++);
               }
            }
            template<class charT, class Traits>
            typename boost::enable_if_c<sizeof(charT) == 4>::type print_unicode_char(std::basic_ostream <charT, Traits>& os, std::uint32_t code_point)
            {
               os.put(code_point);
            }

            template<class charT, class Traits>
            void print_unicode_superscript(std::basic_ostream <charT, Traits>& os, std::basic_string<charT, Traits> digits)
            {
               for (std::size_t i = 0; i < digits.size(); ++i)
               {
                  switch (digits[i])
                  {
                  case '0':
                     print_unicode_char(os, 0x2070);
                     break;
                  case '1':
                     print_unicode_char(os, 0xB9);
                     break;
                  case '2':
                     print_unicode_char(os, 0xB2);
                     break;
                  case '3':
                     print_unicode_char(os, 0xB3);
                     break;
                  case '4':
                     print_unicode_char(os, 0x2074);
                     break;
                  case '5':
                     print_unicode_char(os, 0x2075);
                     break;
                  case '6':
                     print_unicode_char(os, 0x2076);
                     break;
                  case '7':
                     print_unicode_char(os, 0x2076);
                     break;
                  case '8':
                     print_unicode_char(os, 0x2078);
                     break;
                  case '9':
                     print_unicode_char(os, 0x2079);
                     break;
                  case '+':
                     print_unicode_char(os, 0x207A);
                     break;
                  case '-':
                     print_unicode_char(os, 0x207B);
                     break;
                  }
               }
            }

            template<class charT, class Traits>
            void print_unicode_subscript(std::basic_ostream <charT, Traits>& os, std::basic_string<charT, Traits> digits)
            {
               for (std::size_t i = 0; i < digits.size(); ++i)
               {
                  if ((digits[i] >= '0') && (digits[i] <= '9'))
                     print_unicode_char(os, 0x2080 + (digits[i] - '0'));
                  else if (digits[i] == '-')
                     print_unicode_char(os, 0x208B);
                  else if (digits[i] == '+')
                     print_unicode_char(os, 0x208A);
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

         }  //  namespace detail

         enum output_print_t
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

         template <output_print_t Format = text_format, class charT = char, class Traits = std::char_traits<charT> >
         class basic_numeric_printer;
         template <class charT, class Traits>
         class basic_numeric_printer_base;

         template <class charT, class Traits, class Value>
         typename boost::enable_if_c<detail::is_integer_like<Value>::value>::type print(basic_numeric_printer_base<charT, Traits>& os, const Value& value)
         {
            typename basic_numeric_printer_base<charT, Traits>::scoped_prolog s(&os, "integer");
            os.stream() << value;
         }
         template <class charT, class Traits, class Value>
         typename boost::enable_if_c<detail::is_float_like<Value>::value>::type print(basic_numeric_printer_base<charT, Traits>& os, const Value& f)
         {
            using std::isinf;
            using std::isnan;

            typename basic_numeric_printer_base<charT, Traits>::scoped_prolog s(&os, "float");

            if ((isinf)(f))
            {
               if (f < 0)
                  os.stream() << "-";
               os.print_special_character(0x221E);
               return;
            }
            if ((isnan)(f))
            {
               os.print_name("NaN");
               return;
            }
            std::basic_string<charT, Traits> mantissa, exponent;
            os.decompose_float(os.stream(), f, mantissa, exponent);
            os.stream() << mantissa;
            if (exponent.size())
            {
               os.print_times();
               os.stream() << "10";
               os.print_superscript(exponent);
            }
            return;
         }
         template <class charT, class Traits, class Complex>
         typename boost::enable_if_c<detail::is_complex_like<Complex>::value>::type print(basic_numeric_printer_base<charT, Traits>& os, const Complex& f)
         {
            using std::isnan;
            using std::isinf;

            typename basic_numeric_printer_base<charT, Traits>::scoped_prolog s1(&os, "complex");
            typename basic_numeric_printer_base<charT, Traits>::scoped_styling s2(&os);

            if ((isnan)(f.real()) || (isnan)(f.imag()))
            {
               os.print_name("NaN");
               return;
            }
            else if ((isinf)(f.real()) || (isinf)(f.imag()))
            {
               os.print_complex_infinity();
               return;
            }
            else if ((!os.show_zero_components()) && (f.imag() == 0))
            {
               print(os, f.real());
               return;
            }
            else if ((!os.show_zero_components()) && (f.real() == 0))
            {
               print(os, f.imag());
            }
            else
            {
               if (os.parenthesis())
                  os.stream() << "(";
               print(os, f.real());
               typename Complex::value_type i(f.imag());
               bool isneg = i < 0;
               if (isneg)
               {
                  i = -i;
                  os.stream() << " - ";
               }
               else
                  os.stream() << " + ";
               print(os, i);
            }
            os.print_imaginary_unit();
            if (os.parenthesis() && !((!os.show_zero_components()) && (f.real() == 0)))
               os.stream() << ")";
            return;
         }
         template <class charT, class Traits, class Value>
         typename boost::enable_if_c<detail::is_polynomial_like<Value>::value>::type print(basic_numeric_printer_base<charT, Traits>& os, const Value& f)
         {
            typename basic_numeric_printer_base<charT, Traits>::scoped_prolog s1(&os, "polynomial");
            typename basic_numeric_printer_base<charT, Traits>::scoped_parenthesis s2(&os);
            typename basic_numeric_printer_base<charT, Traits>::scoped_styling s3(&os);

            bool have_first = false;
            for (unsigned i = 0; i <= f.degree(); ++i)
            {
               auto coef = f[i];
               if (os.show_zero_components() || (!detail::iszero(coef)))
               {
                  if (have_first)
                  {
                     if (!detail::isneg(coef))
                        os.stream() << " + ";
                     else
                        os.stream() << " - ";
                  }
                  print(os, !detail::isneg(coef) ? coef : -coef);
                  have_first = true;
                  if (i)
                  {
                     os.print_variable('x');
                     if (i > 1)
                     {
                        std::basic_stringstream<charT, Traits> ss;
                        ss << i;
                        std::basic_string<charT, Traits> s = ss.str();
                        os.print_superscript(s);
                     }
                  }
               }
            }
         }
         template <class charT, class Traits, class Value>
         typename boost::enable_if_c<detail::is_rational_like<Value>::value>::type print(basic_numeric_printer_base<charT, Traits>& os, const Value& rat)
         {
            typename basic_numeric_printer_base<charT, Traits>::scoped_prolog s1(&os, "rational");
            typename basic_numeric_printer_base<charT, Traits>::scoped_styling s3(&os);

            auto num = numerator(rat);
            auto denom = denominator(rat);

            if ((num == 0) || (denom == 1))
            {
               print(os, num);
               return;
            }
            else
            {
               if (detail::isneg(num))
               {
                  os.stream() << "-";
                  num = -num;
               }
               if (os.parenthesis())
                  os.stream() << "(";
               {
                  typename basic_numeric_printer_base<charT, Traits>::scoped_parenthesis s2(&os);

                  auto str1 = os.part_as_string(num);
                  auto str2 = os.part_as_string(denom);

                  os.print_fraction(str1, str2);
               }
               if (os.parenthesis())
                  os.stream() << ")";
            }
         }

         template <class charT, class Traits, class Value>
         typename boost::enable_if_c<detail::is_interval_like<Value>::value>::type print(basic_numeric_printer_base<charT, Traits>& os, const Value& val)
         {
            typename basic_numeric_printer_base<charT, Traits>::scoped_prolog s1(&os, "interval");
            typename basic_numeric_printer_base<charT, Traits>::scoped_styling s3(&os);

            if (singleton(val))
            {
               print(os, lower(val));
               return;
            }
            else
            {
               auto med = median(val);
               auto rad = width(val);
               rad /= 2;

               if (os.parenthesis())
                  os.stream() << "(";
               {
                  typename basic_numeric_printer_base<charT, Traits>::scoped_parenthesis s2(&os);
                  print(os, med);
                  os.stream() << " ";
                  os.print_special_character(0xB1);
                  os.stream() << " ";
                  print(os, rad);
               }
               if (os.parenthesis())
                  os.stream() << ")";
            }
         }

         template <class Stream, class Value>
         using printable_t = decltype(print(std::declval<Stream&>(), std::declval<const Value&>()));

         template <class Stream, class Value>
         struct is_printable
         {
            static const bool value = boost::is_detected_v<printable_t, Stream, Value>;
         };


         template <class charT, class Traits>
         class basic_numeric_printer_base
         {
            styling_level_t styling_level;
            imaginary_i_t i_style;
            multiplyer_t  multiply_style;
            std::size_t m_parenthesis;
            bool m_show_zero_components;
            bool m_use_unicode;
            std::basic_ostream<charT, Traits>* p_stream;
         public:
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
            basic_numeric_printer_base(std::basic_ostream<charT, Traits>& os) 
               : styling_level(full_styling), i_style(upright_i), multiply_style(multiply_times), 
                 m_parenthesis(0), m_show_zero_components(false), m_use_unicode(true), p_stream(&os) {}
            basic_numeric_printer_base(std::basic_ostream<charT, Traits>& os, const basic_numeric_printer_base& o) 
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
               basic_numeric_printer_base<charT, Traits>* m_formatter;
               scoped_parenthesis(basic_numeric_printer_base<charT, Traits>* formatter) : m_formatter(formatter)
               {
                  m_formatter->parenthesis(1 + m_formatter->parenthesis());
               }
               ~scoped_parenthesis()
               {
                  m_formatter->parenthesis(m_formatter->parenthesis() - 1);
               }
            };
            struct scoped_styling
            {
               basic_numeric_printer_base<charT, Traits>* m_formatter;
               styling_level_t saved_style;
               scoped_styling(basic_numeric_printer_base<charT, Traits>* formatter) : m_formatter(formatter)
               {
                  saved_style = m_formatter->styling();
                  m_formatter->styling(no_styling);
               }
               ~scoped_styling()
               {
                  m_formatter->styling(saved_style);
               }
            };
            struct scoped_prolog
            {
               basic_numeric_printer_base<charT, Traits>* m_formatter;
               scoped_prolog(basic_numeric_printer_base<charT, Traits>* formatter, const char* number_type) : m_formatter(formatter)
               {
                  m_formatter->print_prolog(number_type);
               }
               ~scoped_prolog()
               {
                  m_formatter->print_epilog();
               }
            };

            virtual void print_prolog(const char*) {}
            virtual void print_epilog() {}
            virtual void print_special_character(boost::uint32_t unicode_value)
            {
               if (use_unicode())
                  detail::print_unicode_char(this->stream(), unicode_value);
               else
               {
                  switch (unicode_value)
                  {
                  case 0x221E:
                     stream() << "INFINITY";
                     break;
                  case 0xD7:
                  case 0x22C5:
                     stream() << "x";
                     break;
                  case 0x2044:
                     stream() << "/";
                  default:
                     throw std::runtime_error("Unsuported Unicode character encountered in plain text output");
                  }
               }
            }
            virtual void print_times()
            {
               if (multiply() == multiply_times)
                  print_special_character(0xD7);
               else if (multiply() == multiply_dot)
                  print_special_character(0x22C5);
               else
                  print_name("x");
            }
            virtual void print_superscript(const std::basic_string<charT, Traits>& s)
            {
               if (use_unicode())
                  detail::print_unicode_superscript(stream(), s);
               else
                  stream() << "^" << s;
            }
            virtual void print_imaginary_unit()
            {
               if ((use_unicode()) && (this->imaginary_style() == doublestruck_i))
                  detail::print_unicode_char(stream(), 0x2148);
               else if ((use_unicode()) && (this->imaginary_style() == slanted_i))
                  detail::print_unicode_char(stream(), 0x1D456);
               else
                  stream().put(stream().widen('i'));
            }
            virtual void print_complex_infinity()
            {
               if (use_unicode())
               {
                  detail::print_unicode_char(stream(), 0x221E);
               }
               else
                  stream() << "COMPLEX INFINITY";
            }
            virtual void print_variable(charT c)
            {
               stream().put(c);
            }
            virtual void print_name(const std::basic_string<charT, Traits>& s)
            {
               stream() << s;
            }
            virtual void print_fraction(const std::basic_string<charT, Traits>& s1, const std::basic_string<charT, Traits>& s2)
            {
               print_superscript(s1);
               print_special_character(0x2044);
               print_subscript(s2);
            }
            virtual void print_subscript(const std::basic_string<charT, Traits>& s)
            {
               detail::print_unicode_subscript(stream(), s);
            }

            bool latex_as_equation()const { return false; }
            void latex_as_equation(bool) {}

            template <class Value>
            std::basic_string<charT, Traits> part_as_string(const Value& value)
            {
               std::basic_ostringstream<charT, Traits> ss;
               ss.copyfmt(stream());
               ss.imbue(stream().getloc());
               ss.width(0);
               basic_numeric_printer_base fmt(ss, *this);
               print(fmt, value);
               return ss.str();
            }
         };

         template <class charT, class Traits>
         class basic_numeric_printer<text_format, charT, Traits> : public basic_numeric_printer_base<charT, Traits>
         {
         public:
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os) : basic_numeric_printer_base<charT, Traits>(os) {}
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os, basic_numeric_printer const& fmt) : basic_numeric_printer_base<charT, Traits>(os, fmt) {}
         };

         template <class charT, class Traits>
         class basic_numeric_printer<docbook_format, charT, Traits> : public basic_numeric_printer_base<charT, Traits>
         {
            friend class basic_numeric_printer_base<charT, Traits>;
            bool need_epilog;
         public:
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os) : basic_numeric_printer_base<charT, Traits>(os), need_epilog(true) {}
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os, basic_numeric_printer const& fmt) : basic_numeric_printer_base<charT, Traits>(os, fmt), need_epilog(true) {}

            virtual void print_prolog(const char* name)
            {
               if (name && *name)
               {
                  if (this->styling() >= minimal_styling)
                     this->stream() << "<phrase role=\"number\">";
                  if (this->styling() > minimal_styling)
                     this->stream() << "<phrase role=\"" << name << "\">";
               }
               else
                  need_epilog = false;
            }
            virtual void print_epilog()
            {
               if (need_epilog)
               {
                  if (this->styling() > minimal_styling)
                     this->stream() << "</phrase>";
                  if (this->styling() >= minimal_styling)
                     this->stream() << "</phrase>";
               }
            }

            virtual void print_special_character(boost::uint32_t unicode_value)
            {
               std::ios_base::fmtflags f = this->stream().flags() & std::ios_base::basefield;
               this->stream() << std::hex << "&#x" << unicode_value << ";";
               this->stream().setf(f, std::ios_base::basefield);
               BOOST_ASSERT(f == (this->stream().flags() & std::ios_base::basefield));
            }
            virtual void print_superscript(const std::basic_string<charT, Traits>& s)
            {
               this->stream() << "<superscript>" << s << "</superscript>";
            }
            virtual void print_imaginary_unit()
            {
               if (this->styling() >= minimal_styling)
                  this->stream() << "<phrase role=\"imaginary_i\">";
               if (this->imaginary_style() == doublestruck_i)
                  this->stream() << "&#x2148;";
               else if (this->imaginary_style() == slanted_i)
                  this->stream() << "<emphasis>i</emphasis>";
               else
                  this->stream().put(this->stream().widen('i'));
               if (this->styling() >= minimal_styling)
                  this->stream() << "</phrase>";
            }
            virtual void print_complex_infinity()
            {
               this->stream() << "&#8734;&#x0303;";
            }
            virtual void print_variable(charT c)
            {
               this->stream() << "<emphasis>";
               this->stream().put(c);
               this->stream() << "</emphasis>";
            }
            virtual void print_subscript(const std::basic_string<charT, Traits>& s)
            {
               this->stream() << "<subscript>" << s << "</subscript>";
            }
         };

         template <class charT, class Traits>
         class basic_numeric_printer<latex_format, charT, Traits> : public basic_numeric_printer_base<charT, Traits>
         {
            friend class basic_numeric_printer_base<charT, Traits>;

            bool is_latex_as_equation;
            unsigned inside_equation;

         public:
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os) : basic_numeric_printer_base<charT, Traits>(os), is_latex_as_equation(true), inside_equation(false) {}
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os, basic_numeric_printer const& fmt) : basic_numeric_printer_base<charT, Traits>(os, fmt), is_latex_as_equation(fmt.is_latex_as_equation), inside_equation(fmt.inside_equation) {}

            void latex_as_equation(bool b) { is_latex_as_equation = b; }
            bool latex_as_equation()const { return is_latex_as_equation; }

            virtual void print_prolog(const char*)
            {
               if (latex_as_equation())
               {
                  if(inside_equation++ == 0)
                     this->stream().put(this->stream().widen('$'));
               }
            }
            virtual void print_epilog()
            {
               if (latex_as_equation())
               {
                  if(--inside_equation == 0)
                     this->stream().put(this->stream().widen('$'));
               }
            }

            virtual void print_special_character(boost::uint32_t unicode_value)
            {
               if (!inside_equation)
                  this->stream() << "$";
               switch (unicode_value)
               {
               case 0x221E:
                  this->stream() << "\\infty";
                  break;
               case 0xD7:
                  this->stream() << "\\times";
                  break;
               case 0x22C5:
                  this->stream() << "\\cdot";
                  break;
               case 0xB1:
                  this->stream() << "\\pm";
                  break;
               case 0x2213:
                  this->stream() << "\\mp";
                  break;
               case 0x221A:
                  this->stream() << "\\surd";
                  break;
               case 0x3B4:
                  this->stream() << "\\delta";
                  break;
               case 0x394:
                  this->stream() << "\\Delta";
                  break;
               case 0x1D70B:
                  this->stream() << "\\pi";
                  break;
               case 0x03C0:
                  this->stream() << "\\pi";
                  break;
               case 0x1D70E:
                  this->stream() << "\\sigma";
                  break;
               case 0x3c3:
                  this->stream() << "\\sigma";
                  break;
               default:
                  throw std::runtime_error("Unsuported Unicode character encountered in LaTex output");
               }
               if (!inside_equation)
                  this->stream() << "$";
            }
            virtual void print_superscript(const std::basic_string<charT, Traits>& s)
            {
               if (inside_equation)
                  this->stream() << "^{" << s << "}";
               else
                  this->stream() << "\\textsuperscript{" << s << "}";
            }
            virtual void print_imaginary_unit()
            {
               if (this->imaginary_style() == doublestruck_i)
                  BOOST_THROW_EXCEPTION(std::runtime_error("Doublestruck imaginary i not implemented for LaTex output."));
               else if ((this->imaginary_style() == upright_i) && inside_equation)
                  this->stream() << "\\mathrm{i}";
               else if ((this->imaginary_style() == slanted_i) && !inside_equation)
                  this->stream() << "\\textit{i}";
               else
                  this->stream().put(this->stream().widen('i'));
            }
            virtual void print_complex_infinity()
            {
               if (!inside_equation)
                  this->stream().put(this->stream().widen('$'));
               this->stream() << "\\tilde{\\infty}";
               if (!inside_equation)
                  this->stream().put(this->stream().widen('$'));
            }
            virtual void print_variable(charT c)
            {
               if (inside_equation)
                  this->stream().put(c);
               else
               {
                  this->stream() << "\\textit{";
                  this->stream().put(c);
                  this->stream() << "}";
               }
            }
            virtual void print_name(const std::basic_string<charT, Traits>& s)
            {
               if (!inside_equation)
                  this->stream() << s;
               else
               {
                  this->stream() << "\\mathrm{";
                  this->stream() << s;
                  this->stream() << "}";
               }
            }
            virtual void print_fraction(const std::basic_string<charT, Traits>& s1, const std::basic_string<charT, Traits>& s2)
            {
               if (inside_equation)
               {
                  this->stream() << "\\frac{" << s1 << "}{" << s2 << "}";
               }
               else
               {
                  this->print_superscript(s1);
                  this->stream() << "/";
                  this->print_subscript(s2);
               }
            }
            virtual void print_subscript(const std::basic_string<charT, Traits>& s)
            {
               if (inside_equation)
                  this->stream() << "_{" << s << "}";
               else
                  this->stream() << "\\textsubscript{" << s << "}";
            }

         };

         template <class charT, class Traits>
         class basic_numeric_printer<html_format, charT, Traits> : public basic_numeric_printer_base<charT, Traits>
         {
            friend class basic_numeric_printer_base<charT, Traits>;
            bool need_epilog;
         public:
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os) : basic_numeric_printer_base<charT, Traits>(os), need_epilog(true){}
            basic_numeric_printer(std::basic_ostream<charT, Traits>& os, basic_numeric_printer const & fmt) : basic_numeric_printer_base<charT, Traits>(os, fmt), need_epilog(true){}

            virtual void print_prolog(const char* name)
            {
               if (name && *name)
               {
                  if (this->styling() >= minimal_styling)
                     this->stream() << "<span class=\"number\">";
                  if (this->styling() > minimal_styling)
                     this->stream() << "<span class=\"" << name << "\">";
               }
               else
                  need_epilog = false;
            }
            virtual void print_epilog()
            {
               if (need_epilog)
               {
                  if (this->styling() > minimal_styling)
                     this->stream() << "</span>";
                  if (this->styling() >= minimal_styling)
                     this->stream() << "</span>";
               }
            }

            virtual void print_special_character(boost::uint32_t unicode_value)
            {
               std::ios_base::fmtflags f = this->stream().flags() & std::ios_base::basefield;
               this->stream() << std::hex << "&#x" << unicode_value << ";";
               this->stream().setf(f, std::ios_base::basefield);
               BOOST_ASSERT(f == (this->stream().flags() & std::ios_base::basefield));
            }
            virtual void print_superscript(const std::basic_string<charT, Traits>& s)
            {
               this->stream() << "<sup>" << s << "</sup>";
            }
            virtual void print_imaginary_unit()
            {
               if (this->styling() >= minimal_styling)
                  this->stream() << "<span class=\"imaginary_i\">";
               if (this->imaginary_style() == doublestruck_i)
                  this->stream() << "&#x2148;";
               else if (this->imaginary_style() == slanted_i)
                  this->stream() << "<I>i</I>";
               else
                  this->stream().put(this->stream().widen('i'));
               if (this->styling() >= minimal_styling)
                  this->stream() << "</span>";
            }
            virtual void print_complex_infinity()
            {
               this->stream() << "&#8734;&#x0303;";
            }
            virtual void print_variable(charT c)
            {
               this->stream() << "<I>";
               this->stream().put(c);
               this->stream() << "</I>";
            }
            virtual void print_subscript(const std::basic_string<charT, Traits>& s)
            {
               this->stream() << "<sub>" << s << "</sub>";
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

         template <output_print_t Format, class charT, class Traits, class Value>
         typename boost::enable_if_c<is_printable<basic_numeric_printer<Format, charT, Traits>, Value>::value, basic_numeric_printer<Format, charT, Traits>&>::type operator << (basic_numeric_printer<Format, charT, Traits>& os, const Value& i)
         {
            std::size_t w = (std::size_t)os.stream().width();
            if (w)
            {
               std::basic_string<charT, Traits> s = os.part_as_string(i);
               std::size_t len = detail::unicode_character_len(s.c_str());
               if (len < w)
               {
                  auto pos = os.stream().flags() & std::ios_base::adjustfield;
                  s.insert(pos == std::ios_base::left ? s.end() : s.begin(), w - len, os.stream().fill());
               }
               os.stream() << s;
            }
            else
               print(os, i);
            return os;
         }

         template <output_print_t Format, class charT, class Traits, class T>
         typename boost::disable_if_c<is_printable<basic_numeric_printer<Format, charT, Traits>, T>::value, basic_numeric_printer<Format, charT, Traits>&>::type operator<< (basic_numeric_printer<Format, charT, Traits>& os, const T& unhandled)
         {
            os.stream() << unhandled;
            return os;
         }
         //
         // Manipulators:
         //
         template <output_print_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, styling_level_t s)
         {
            os.styling(s);
            return os;
         }

         template <output_print_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, multiplyer_t s)
         {
            os.multiply(s);
            return os;
         }

         template <output_print_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, latex_t l)
         {
            os.latex_as_equation(l == latex_as_equation);
            return os;
         }

         template <output_print_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, const imaginary_i_t i)
         {
            os.imaginary_style(i);
            return os;
         }

         template <output_print_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, const zero_component_t z)
         {
            os.show_zero_components(z == ::boost::math::tools::show_zero_components ? true : false);
            return os;
         }

         template <output_print_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator << (basic_numeric_printer<Format, charT, Traits>& os, const unicode_text_t z)
         {
            os.use_unicode(z == ::boost::math::tools::unicode_text_output ? true : false);
            return os;
         }

         //
         // Output formatters:
         //
         template <output_print_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator<< (basic_numeric_printer<Format, charT, Traits>& os, std::basic_ostream<charT, Traits>& (*pf)(std::basic_ostream<charT, Traits>&))
         {
            pf(os.stream());
            return os;
         }
         template <output_print_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator<< (basic_numeric_printer<Format, charT, Traits>& os, std::basic_ios<charT, Traits>& (*pf)(std::basic_ios<charT, Traits>&))
         {
            pf(os.stream());
            return os;
         }
         template <output_print_t Format, class charT, class Traits>
         basic_numeric_printer<Format, charT, Traits>& operator<< (basic_numeric_printer<Format, charT, Traits>& os, std::ios_base& (*pf)(std::ios_base&))
         {
            pf(os.stream());
            return os;
         }
      }
   }
}

#endif // BOOST_MATH_TOOLS_FORMATTING_HPP

