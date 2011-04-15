/** nonfinite_num_facet_serialization.cpp
 *
 * Copyright (c) 2011 Francois Mauger
 * Copyright (c) 2011 Paul A. Bristow
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt
 * or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * This    sample    program    illustrates    how    to    use    the
 * `boost/math/nonfinite_num_facets.hpp'  material  from the  original
 * Floating Point  Utilities contribution by  Johan Rade.  Here  it is
 * shown  how  non  finite  floating  number  can  be  serialized  and
 * deserialized from  I/O streams and/or Boost  text/XML archives.  It
 * produces two archives stored in `test.txt' and `test.xml' files.
 *
 * Tested with Boost 1.44, gcc 4.4.1, Linux/i686 (32bits)
 *
 */


#ifdef _MSC_VER
#   pragma warning(push)

#   pragma warning(disable : 4224) // formal parameter 'version' was previously defined as a type
#   pragma warning(disable : 4100) // unreferenced formal parameter
#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>

#include <boost/cstdint.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/codecvt_null.hpp>

// from Johan Rade Floating Point Utilities :
#include <boost/math/special_functions/nonfinite_num_facets.hpp>

/* A class with a float and a double */
struct foo
{
  foo () : fvalue (3.141593F), dvalue (3.14159265358979) {}
  // set the values at -infinity :
  void minus_infinity ()
  {
    fvalue = -std::numeric_limits<float>::infinity ();
    dvalue = -std::numeric_limits<double>::infinity ();
    return;
  }
  // set the values at +infinity :
  void plus_infinity ()
  {
    fvalue = +std::numeric_limits<float>::infinity ();
    dvalue = +std::numeric_limits<double>::infinity ();
    return;
  }
  // set the values at NaN :
  void nan ()
  {
    fvalue = +std::numeric_limits<float>::quiet_NaN ();
    dvalue = +std::numeric_limits<double>::quiet_NaN ();
    return;
  }
  // print :
  void print (std::ostream & a_out, const std::string & a_title)
  {
    if (a_title.empty ()) a_out << "foo";
    else a_out << a_title;
    a_out << " : " << std::endl;
    a_out << "|-- " << "fvalue = ";
    a_out.precision (7);
    a_out << fvalue << std::endl;
    a_out << "`-- " << "dvalue = ";
    a_out.precision (15);
    a_out << dvalue << std::endl;
    return;
  }

  // I/O operators :
  friend std::ostream & operator<< (std::ostream & a_out, const foo & a_foo);
  friend std::istream & operator>> (std::istream & a_in, foo & a_foo);

  // Boost serialization :
  template <class Archive>
  void serialize (Archive & ar, int version)
  {
    ar & BOOST_SERIALIZATION_NVP (fvalue);
    ar & BOOST_SERIALIZATION_NVP (dvalue);
    return;
  }

  // Attributes :
  float  fvalue; // single precision floating number
  double dvalue; // double precision floating number
};

std::ostream & operator<< (std::ostream & a_out, const foo & a_foo)
{
  a_out << "(" << a_foo.fvalue << ";"   << a_foo.dvalue << ")";
  return a_out;
}

std::istream & operator>> (std::istream & a_in, foo & a_foo)
{
  char c = 0;
  a_in.get (c);
  if (c != '(')
  {
    std::cerr << "ERROR: operator>> No ( " << std::endl;
    a_in.setstate(std::ios::failbit);
    return a_in;
  }
  float f;
  a_in >> std::ws >> f;
  if (! a_in)
  {
    return a_in;
  }
  a_in >> std::ws;
  a_in.get (c);
  if (c != ';')
  {
    std::cerr << "ERROR: operator>> c='" << c << "'" << std::endl;
    std::cerr << "ERROR: operator>> No ; " << std::endl;
    a_in.setstate(std::ios::failbit);
    return a_in;
  }
  double d;
  a_in >> std::ws >> d;
  if (! a_in)
  {
    return a_in;
  }
  a_in >> std::ws;
  a_in.get (c);
  if (c != ')')
  {
    std::cerr << "ERROR: operator>> No ) " << std::endl;
    a_in.setstate(std::ios::failbit);
    return a_in;
  }
  a_foo.fvalue = f;
  a_foo.dvalue = d;
  return a_in;
}

int main (void)
{
  std::clog << "Hello Booster !" << std::endl
      << "This is the `test_nonfinite_num_facets_1.cpp' sample program !" << std::endl;

  std::locale the_default_locale (std::locale::classic (),
          new boost::archive::codecvt_null<char>);

  {
    std::clog << "Write to some string buffer..." << std::endl;
    foo f0;
    foo f1; f1.minus_infinity ();
    foo f2; f2.plus_infinity ();
    foo f3; f3.nan ();

    f0.print (std::clog, "f0");
    f1.print (std::clog, "f1");
    f2.print (std::clog, "f2");
    f3.print (std::clog, "f3");

    std::ostringstream oss;
    std::locale the_out_locale (the_default_locale, new boost::math::nonfinite_num_put<char>);
    oss.imbue (the_out_locale);
    oss.precision (15);
    oss << f0 << f1 << f2 << f3;
    std::clog << "Output is: `" << oss.str () << "'" << std::endl;
    std::clog << "Done." << std::endl;
  }

  {
    std::clog << "Read from to some string buffer..." << std::endl;
    std::string the_string = "(3.14159;3.14159)(-inf;-inf)(inf;inf)(nan;nan)";
    std::clog << "Input is: `" << the_string << "'" << std::endl;

    std::locale the_in_locale (the_default_locale, new boost::math::nonfinite_num_get<char>);
    std::istringstream iss (the_string);
    iss.imbue (the_in_locale);

    foo f0, f1, f2, f3;
    iss >> f0 >> f1 >> f2 >> f3;
    if (! iss)
    {
      std::cerr << "Format error !" << std::endl;
    }
    else
    {
      std::cerr << "Success !" << std::endl;
      f0.print (std::clog, "f0");
      f1.print (std::clog, "f1");
      f2.print (std::clog, "f2");
      f3.print (std::clog, "f3");
    }
    std::clog << "Done." << std::endl;
  }

  {
    std::clog << "Serialize (boost text archive)..." << std::endl;
    foo f0;
    foo f1; f1.minus_infinity ();
    foo f2; f2.plus_infinity ();
    foo f3; f3.nan ();

    f0.print (std::clog, "f0");
    f1.print (std::clog, "f1");
    f2.print (std::clog, "f2");
    f3.print (std::clog, "f3");

    std::locale the_out_locale (the_default_locale, new boost::math::nonfinite_num_put<char>);
    std::ofstream fout ("test.txt");
    fout.imbue (the_out_locale);
    boost::archive::text_oarchive toar (fout, boost::archive::no_codecvt);

    toar & f0;
    toar & f1;
    toar & f2;
    toar & f3;
    std::clog << "Done." << std::endl;
  }

  {
    std::clog << "Deserialize (boost text archive)..." << std::endl;
    std::locale the_in_locale (the_default_locale, new boost::math::nonfinite_num_get<char>);
    std::ifstream fin ("test.txt");
    fin.imbue (the_in_locale);
    boost::archive::text_iarchive tiar (fin, boost::archive::no_codecvt);
    foo f0, f1, f2, f3;

    tiar & f0;
    tiar & f1;
    tiar & f2;
    tiar & f3;

    f0.print (std::clog, "f0");
    f1.print (std::clog, "f1");
    f2.print (std::clog, "f2");
    f3.print (std::clog, "f3");

    std::clog << "Done." << std::endl;
  }

  {
    std::clog << "Serialize (boost XML archive)..." << std::endl;
    foo f0;
    foo f1; f1.minus_infinity ();
    foo f2; f2.plus_infinity ();
    foo f3; f3.nan ();

    f0.print (std::clog, "f0");
    f1.print (std::clog, "f1");
    f2.print (std::clog, "f2");
    f3.print (std::clog, "f3");

    std::locale the_out_locale (the_default_locale, new boost::math::nonfinite_num_put<char>);
    std::ofstream fout ("test.xml");
    fout.imbue (the_out_locale);
    boost::archive::xml_oarchive xoar (fout, boost::archive::no_codecvt);

    xoar & BOOST_SERIALIZATION_NVP (f0);
    xoar & BOOST_SERIALIZATION_NVP (f1);
    xoar & BOOST_SERIALIZATION_NVP (f2);
    xoar & BOOST_SERIALIZATION_NVP (f3);
    std::clog << "Done." << std::endl;
  }

  {
    std::clog << "Deserialize (boost XML archive)..." << std::endl;
    std::locale the_in_locale (the_default_locale, new boost::math::nonfinite_num_get<char>);
    std::ifstream fin ("test.xml");
    fin.imbue (the_in_locale);
    boost::archive::xml_iarchive xiar (fin, boost::archive::no_codecvt);
    foo f0, f1, f2, f3;

    xiar & BOOST_SERIALIZATION_NVP (f0);
    xiar & BOOST_SERIALIZATION_NVP (f1);
    xiar & BOOST_SERIALIZATION_NVP (f2);
    xiar & BOOST_SERIALIZATION_NVP (f3);

    f0.print (std::clog, "f0");
    f1.print (std::clog, "f1");
    f2.print (std::clog, "f2");
    f3.print (std::clog, "f3");

    std::clog << "Done." << std::endl;
  }

  std::clog << "Bye ! " << std::endl;
  return 0;
}

 /* end of test_nonfinite_num_facets_1.cpp */
