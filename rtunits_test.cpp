/*
 * Copyright (c) 2019, Ben Barsdell. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * * Neither the name of the copyright holder nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// TODO: Test all methods and operators of Dimensions and Quantity.

#define RTUNITS_USE_EXCEPTIONS 1
// TODO: Makes the CI run tests twice, once with this on and once off.
//#define RTUNITS_USE_EXCEPTIONS 0
#include "rtunits.hpp"

#include <unordered_set>

#include "gtest/gtest.h"

using namespace rtunits;

constexpr const Dimensions::BaseDimension kBaseDimensionCount =
    Dimensions::kBaseDimensionCount;

TEST(UnitsTest, DimensionsArithmetic) {
  const auto none = Dimensions::None();
  EXPECT_FALSE(static_cast<bool>(none));

  for (int d = 0; d < kBaseDimensionCount; ++d) {
    Dimensions dim(static_cast<Dimensions::BaseDimension>(d));
    EXPECT_TRUE(static_cast<bool>(dim));
    EXPECT_EQ(dim / dim, none);
    EXPECT_EQ(dim.squared() / dim.squared(), none);
    EXPECT_EQ(dim.cubed() / dim.cubed(), none);
    EXPECT_EQ(sqrt(dim.squared()), dim);
    EXPECT_EQ(cbrt(dim.cubed()), dim);
    EXPECT_EQ(dim.power(4).fractional_power(4), dim);
  }
}

#define EXPECT_QUANT_EQ(_lhs, _rhs)                       \
  do {                                                    \
    const auto& lhs = (_lhs);                             \
    const auto& rhs = (_rhs);                             \
    EXPECT_EQ((lhs).dimensions(), (rhs).dimensions());    \
    EXPECT_DOUBLE_EQ((lhs).value_si(), (rhs).value_si()); \
  } while (0)

TEST(UnitsTest, QuantityArithmetic) {
  const auto unity = Quantity64(1.);
  for (const auto& unit_item : Quantity64::unit_symbol_map()) {
    const Quantity64& q = unit_item.second;
    EXPECT_QUANT_EQ(q * -1., -q);
    EXPECT_QUANT_EQ(q * +1., +q);
    EXPECT_QUANT_EQ(q / q, unity);
    EXPECT_QUANT_EQ(q.squared() / q.squared(), unity);
    EXPECT_QUANT_EQ(q.cubed() / q.cubed(), unity);
    EXPECT_QUANT_EQ(sqrt(q.squared()), q);
    EXPECT_QUANT_EQ(cbrt(q.cubed()), q);
    EXPECT_QUANT_EQ(q.power(4).fractional_power(4), q);
    EXPECT_QUANT_EQ(fmod(q * 1.25, q), q * 0.25);
    Quantity64 qmut = q;
    qmut *= q;
    qmut += q * q;
    qmut -= q * q;
    qmut /= q;
    qmut *= 2;
    EXPECT_QUANT_EQ(qmut, q * 2.);
  }
}

TEST(UnitsTest, QuantityComparisons) {
  for (const auto& unit_item : Quantity64::unit_symbol_map()) {
    const Quantity64& q = unit_item.second;
    EXPECT_TRUE(q == q);
    EXPECT_FALSE(q != q);
    EXPECT_FALSE(q < q);
    EXPECT_TRUE(q <= q);
    EXPECT_FALSE(q > q);
    EXPECT_TRUE(q >= q);
  }
}

TEST(UnitsTest, DimensionsComparisons) {
  const auto d = Dimensions::Length();
  EXPECT_TRUE(d == d);
  EXPECT_FALSE(d != d);
  EXPECT_FALSE(d < d);
  EXPECT_TRUE(d <= d);
  EXPECT_FALSE(d > d);
  EXPECT_TRUE(d >= d);
}

TEST(UnitsTest, DimensionsToString) {
  const auto none = Dimensions::None();
  const auto length = Dimensions::Length();
  const auto mass = Dimensions::Mass();
  const auto time = Dimensions::Time();
  const auto current = Dimensions::Current();
  const auto temperature = Dimensions::Temperature();
  const auto luminosity = Dimensions::Luminosity();
  const auto amount = Dimensions::Amount();
  EXPECT_TRUE(none.to_string().empty());
  EXPECT_EQ(none.to_string(true), "<none>");
  EXPECT_EQ(length.to_string(), "L");
  EXPECT_EQ(mass.to_string(), "M");
  EXPECT_EQ(time.to_string(), "T");
  EXPECT_EQ(current.to_string(), "I");
  EXPECT_EQ(temperature.to_string(), "\u03F4" /*capital Theta*/);
  EXPECT_EQ(luminosity.to_string(), "J");
  EXPECT_EQ(amount.to_string(), "N");

  const auto time_squared = time.squared();
  EXPECT_EQ(time_squared.to_string(), "T^2");
  const auto inverse_time = time.reciprocal();
  EXPECT_EQ(inverse_time.to_string(), "T^-1");

  const auto amount_per_time_squared = amount / time.squared();
  EXPECT_EQ(amount_per_time_squared.to_string(), "T^-2 N");

  Dimensions::SymbolArray si_symbols = {"m", "kg", "s", "A", "K", "cd", "mol"};
  EXPECT_EQ(amount_per_time_squared.to_string(si_symbols), "s^-2 mol");
}

TEST(UnitsTest, Hashable) {
  std::unordered_set<Dimensions> dset;
  dset.insert(Dimensions::Length());
  dset.insert(Dimensions::Current());
  dset.insert(Dimensions::Time() * Dimensions::Mass());
  EXPECT_EQ(dset.size(), 3);

  std::unordered_set<Quantity64> qset;
  qset.insert(Quantity64(100, "MHz"));
  qset.insert(Quantity64(10, "s"));
  qset.insert(Quantity64(1, "uF"));
  EXPECT_EQ(qset.size(), 3);
}

TEST(UnitsTest, Sortable) {
  std::set<Dimensions> dset;
  dset.insert(Dimensions::Length());
  dset.insert(Dimensions::Current());
  dset.insert(Dimensions::Time() * Dimensions::Mass());
  EXPECT_EQ(dset.size(), 3);

  std::set<Quantity64> qset;
  qset.insert(Quantity64::second());
  qset.insert(Quantity64::minute());
  qset.insert(Quantity64::hour());
  EXPECT_EQ(qset.size(), 3);

  // Can't mix quantities of different dimensions in an ordered set.
#if RTUNITS_USE_EXCEPTIONS
  EXPECT_THROW(qset.insert(Quantity64::meter()), QuantityError);
#endif
}

#define EXPECT_PARSE_Q(symbol, quantity)    \
  EXPECT_EQ(Quantity64(symbol), quantity)   \
      << "where " #symbol " = \"" << symbol \
      << "\", error = " << fabs(Quantity64(symbol) - quantity);

#if RTUNITS_USE_EXCEPTIONS
#define EXPECT_PARSE_THROW_Q(symbol)                         \
  EXPECT_THROW((void)Quantity64(symbol), QuantityParseError) \
      << "where " #symbol " = \"" << symbol << "\""
#else
#define EXPECT_PARSE_THROW_Q(symbol)
#endif

void test_basic_parsing_with_div(const std::string& symbol,
                                 const Quantity64& quantity) {
  EXPECT_PARSE_Q(symbol, quantity);
  for (const auto& dsym : {"/", " / "}) {
    EXPECT_PARSE_Q(dsym + symbol, quantity.reciprocal());
  }
}

void test_basic_parsing(const std::string& symbol, const Quantity64& quantity) {
  EXPECT_PARSE_Q(symbol, quantity);
  for (const auto& esym : {"", "^", "**", " ", " ^ ", " ** "}) {
    test_basic_parsing_with_div(symbol + esym + "0", 1.);
    test_basic_parsing_with_div(symbol + esym + "1", quantity);
    test_basic_parsing_with_div(symbol + esym + "2", quantity.squared());

    test_basic_parsing_with_div(symbol + esym + "3", quantity.cubed());
    test_basic_parsing_with_div(symbol + esym + "-1", quantity.reciprocal());
    test_basic_parsing_with_div(symbol + esym + "-2", quantity.power(-2));
    test_basic_parsing_with_div(symbol + esym + "-3", quantity.power(-3));
  }
}

TEST(UnitsTest, BasicParsing) {
  EXPECT_EQ(Quantity64(""), Dimensions::None());

  test_basic_parsing("m", Quantity64::meter());
  test_basic_parsing("s", Quantity64::second());
  test_basic_parsing("A", Quantity64::ampere());
  test_basic_parsing("K", Quantity64::kelvin());
  test_basic_parsing("cd", Quantity64::candela());
  test_basic_parsing("mol", Quantity64::mole());
  test_basic_parsing("Hz", Quantity64::hertz());
  test_basic_parsing("rad", Quantity64::radian());
  test_basic_parsing("sr", Quantity64::steradian());
  test_basic_parsing("N", Quantity64::newton());
  test_basic_parsing("Pa", Quantity64::pascal());
  test_basic_parsing("J", Quantity64::joule());
  test_basic_parsing("W", Quantity64::watt());
  test_basic_parsing("C", Quantity64::coulomb());
  test_basic_parsing("V", Quantity64::volt());
  test_basic_parsing("F", Quantity64::farad());
  test_basic_parsing("ohm", Quantity64::ohm());
  test_basic_parsing("S", Quantity64::siemens());
  test_basic_parsing("Wb", Quantity64::weber());
  test_basic_parsing("T", Quantity64::tesla());
  test_basic_parsing("H", Quantity64::henry());
  test_basic_parsing("degC", Quantity64::degree_celsius());
  test_basic_parsing("lm", Quantity64::lumen());
  test_basic_parsing("lx", Quantity64::lux());
  test_basic_parsing("Bq", Quantity64::becquerel());
  test_basic_parsing("Gy", Quantity64::gray());
  test_basic_parsing("Sv", Quantity64::sievert());
  test_basic_parsing("kat", Quantity64::katal());
}

TEST(UnitsTest, UnicodeParsing) {
  test_basic_parsing("\u00B0", Quantity64::degree());
  test_basic_parsing("\u2032", Quantity64::arc_minute());
  test_basic_parsing("\u2033", Quantity64::arc_second());
  test_basic_parsing(
      "\u00B5"
      "as",
      Quantity64::micro() * Quantity64::arc_second());
  test_basic_parsing(
      "\u00B5"
      "\u2033",
      Quantity64::micro() * Quantity64::arc_second());
  test_basic_parsing("\u00C5", Quantity64::angstrom());
  test_basic_parsing("M\u2609", Quantity64::solar_mass());
  test_basic_parsing("R\u2295", Quantity64::earth_radius());

  EXPECT_PARSE_Q(
      "\u00B0"
      "/\u2103",
      Quantity64::degree() / Quantity64::degree_celsius());
}

TEST(UnitsTest, InvalidParsing) {
  for (const auto& unit_item : Quantity64::unit_symbol_map()) {
    const std::string& unit_str = unit_item.first;
    EXPECT_PARSE_THROW_Q(unit_str + unit_str + unit_str);
    EXPECT_PARSE_THROW_Q(unit_str + "xxxxxx");
    EXPECT_PARSE_THROW_Q(unit_str + " xxxxxx");
    EXPECT_PARSE_THROW_Q("xxxxxx" + unit_str);
    EXPECT_PARSE_THROW_Q("xxxxxx " + unit_str);
  }
  EXPECT_PARSE_THROW_Q("kgm");
  EXPECT_PARSE_THROW_Q("'m");
  EXPECT_PARSE_THROW_Q("\"m");
  EXPECT_PARSE_THROW_Q("m ''");
  EXPECT_PARSE_THROW_Q("m \"\"");
  EXPECT_PARSE_THROW_Q("M\u2609M");
}

void test_multi_parsing(std::initializer_list<std::string> symbols,
                        std::initializer_list<Quantity64> quantities,
                        bool check_with_bad_separators = true) {
  const char* middot = "\xB7";
  const char* times = "\xD7";
  for (const auto& leading_ws : {"", "   ", "\t"}) {
    for (const auto& trailing_ws : {"", "   ", "\t"}) {
      for (const auto& sep :
           {" ", "  ", "\t", "*", ",", " * ", " , ", middot, times}) {
        std::string combined_symbol = leading_ws;
        for (const std::string& symbol : symbols) {
          combined_symbol += symbol + sep;
        }
        combined_symbol += trailing_ws;
        Quantity64 combined_quantity(1.);
        for (const Quantity64& quantity : quantities) {
          combined_quantity *= quantity;
        }
        EXPECT_PARSE_Q(combined_symbol, combined_quantity);
      }
      if (check_with_bad_separators) {
        for (const auto& bad_sep : {"", ".", "-", "&", "|"}) {
          std::string combined_symbol = leading_ws;
          for (const std::string& symbol : symbols) {
            combined_symbol += symbol + bad_sep;
          }
          combined_symbol += trailing_ws;
          EXPECT_PARSE_THROW_Q(combined_symbol);
        }
      }
    }
  }
}

TEST(UnitsTest, MultiUnitParsing) {
  test_multi_parsing({"m", "/ s", "/s"},
                     {Quantity64::meter(), Quantity64::second().reciprocal(),
                      Quantity64::second().reciprocal()},
                     false);
  test_multi_parsing({"kg", "m", "s-2"},
                     {Quantity64::kilogram(), Quantity64::meter(),
                      Quantity64::second().squared().reciprocal()});
  auto units = {"m", "s", "P"};
  for (auto unit0 : units) {
    for (auto unit1 : units) {
      for (auto unit2 : units) {
        Quantity64 quant0 = Quantity64::unit_symbol_map().at(unit0);
        Quantity64 quant1 = Quantity64::unit_symbol_map().at(unit1);
        Quantity64 quant2 = Quantity64::unit_symbol_map().at(unit2);
        test_multi_parsing({unit0, unit1, unit2}, {quant0, quant1, quant2});
      }
    }
  }
}

TEST(UnitsTest, PrefixParsing) {
  for (const auto& prefix_item : Quantity64::prefix_symbol_map()) {
    const std::string& prefix_str = prefix_item.first;
    const Quantity64& prefix = prefix_item.second;
    EXPECT_FALSE(static_cast<bool>(prefix.dimensions()));
    for (const auto& unit_item : Quantity64::unit_symbol_map()) {
      const std::string& unit_str = unit_item.first;
      const Quantity64& unit = unit_item.second;
      EXPECT_PARSE_Q(unit_str, unit);
      EXPECT_PARSE_Q(prefix_str + unit_str, prefix * unit);
    }
  }
}
