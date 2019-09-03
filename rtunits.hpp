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

// RT-Units - A single-header library for working with physical quantities
// at runtime. Includes support for parsing and serialization.

// TODO: Consider support for ratios as dimensions using 1 extra bit per dim.
//         E.g., sr == m**2 / m**2 => dim[m] = 2, is_ratio[m] = true.
//         cos(m/m) -> m/m
// TODO: Consider a C API with functions like:
//         int32_t units_power(const char* units, int8_t power,
//                             uint8_t result_size, char* result);
//         int32_t units_multiply(const char* lhs_units, const char* rhs_units,
//                                uint8_t result_size, char* result);
// TODO: Try to make Quantity() constexpr.

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// This macro and the associated helper class come from here:
// https://gist.github.com/oliora/928424f7675d58fadf49c70fdba70d2f
// When evaluated at compile time emits a compilation error if condition is not
// true. Invokes the standard assert at run time.
#define CONSTEXPR_ASSERT(cond)                                                 \
  ((void)((cond) ? 0                                                           \
                 : (detail::constexpr_assert_failed([]() { assert(!#cond); }), \
                    0)))

#ifndef RTUNITS_USE_EXCEPTIONS
#define RTUNITS_USE_EXCEPTIONS 1
#endif

#if RTUNITS_USE_EXCEPTIONS
#include <stdexcept>
#define UNITS_ASSERT(pred, exception) ((pred) ? (void)0 : (throw exception))
#define UNITS_THROW(exception) (throw exception)
#else
#define UNITS_ASSERT(pred, exception) CONSTEXPR_ASSERT(pred)
#define UNITS_THROW(exception) ((void)0)
#endif  // RTUNITS_USE_EXCEPTIONS

namespace rtunits {

namespace detail {

template <class Assert>
inline void constexpr_assert_failed(Assert&& a) noexcept {
  std::forward<Assert>(a)();
}

}  // end namespace detail

#if RTUNITS_USE_EXCEPTIONS
struct DimensionError : public std::runtime_error {
  DimensionError(const std::string& text) : std::runtime_error(text) {}
};

struct QuantityError : public std::runtime_error {
  QuantityError(const std::string& text) : std::runtime_error(text) {}
};

struct QuantityParseError : public QuantityError {
  QuantityParseError(const std::string& text) : QuantityError(text) {}
};
#endif

class Dimensions {
 public:
  enum BaseDimension : uint8_t {
    kLength,
    kMass,
    kTime,
    kCurrent,
    kTemperature,
    kLuminosity,
    kAmount,
    kBaseDimensionCount,
  };

  static constexpr Dimensions None() { return Dimensions(); }
  static constexpr Dimensions Length() { return Dimensions(1); }
  static constexpr Dimensions Mass() { return Dimensions(0, 1); }
  static constexpr Dimensions Time() { return Dimensions(0, 0, 1); }
  static constexpr Dimensions Current() { return Dimensions(0, 0, 0, 1); }
  static constexpr Dimensions Temperature() {
    return Dimensions(0, 0, 0, 0, 1);
  }
  static constexpr Dimensions Luminosity() {
    return Dimensions(0, 0, 0, 0, 0, 1);
  }
  static constexpr Dimensions Amount() {
    return Dimensions(0, 0, 0, 0, 0, 0, 1);
  }

  typedef int8_t value_type;
  typedef std::array<value_type, kBaseDimensionCount> array_type;
  typedef std::tuple<value_type, value_type, value_type, value_type, value_type,
                     value_type, value_type>
      tuple_type;

  constexpr Dimensions() : exponents_() {}

  explicit Dimensions(BaseDimension dim) : Dimensions() { exponents_[dim] = 1; }

  constexpr const array_type& exponents() const { return exponents_; }

  // These constexpr methods technically require C++14, but GCC 5.3 compiles
  // them under C++11 mode without issue.
  constexpr bool operator==(const Dimensions& rhs) const {
#if __cplusplus > 201703L
    return exponents_ == rhs.exponents_;
#endif
    return as_tuple() == rhs.as_tuple();
  }

  constexpr bool operator!=(const Dimensions& rhs) const {
#if __cplusplus > 201703L
    return exponents_ != rhs.exponents_;
#endif
    return as_tuple() != rhs.as_tuple();
  }

  constexpr bool operator<(const Dimensions& rhs) const {
#if __cplusplus > 201703L
    return exponents_ < rhs.exponents_;
#endif
    return as_tuple() < rhs.as_tuple();
  }

  constexpr bool operator<=(const Dimensions& rhs) const {
    return *this == rhs || *this < rhs;
  }

  constexpr bool operator>(const Dimensions& rhs) const {
    return !(*this <= rhs);
  }

  constexpr bool operator>=(const Dimensions& rhs) const {
    return !(*this < rhs);
  }

  // Test if this has any dimensions.
  explicit constexpr operator bool() const {
    return length() || mass() || time() || current() || temperature() ||
           luminosity() || amount();
  }

  // Test if this is a base unit (has exactly 1 unit set).
  constexpr bool is_base() const {
    return bool(length()) + bool(mass()) + bool(time()) + bool(current()) +
               bool(temperature()) + bool(luminosity()) + bool(amount()) ==
           1;
  }

  // TODO: Make constexpr somehow (also outside operators).
  Dimensions& operator*=(const Dimensions& rhs) {
    std::transform(exponents_.begin(), exponents_.end(), rhs.exponents_.begin(),
                   exponents_.begin(), std::plus<value_type>());
    return *this;
  }

  Dimensions& operator/=(const Dimensions& rhs) {
    std::transform(exponents_.begin(), exponents_.end(), rhs.exponents_.begin(),
                   exponents_.begin(), std::minus<value_type>());
    return *this;
  }

  constexpr Dimensions power(value_type n) const {
    return Dimensions(length() * n, mass() * n, time() * n, current() * n,
                      temperature() * n, luminosity() * n, amount() * n);
  }

  constexpr Dimensions reciprocal() const { return power(-1); }

  constexpr Dimensions squared() const { return power(2); }

  constexpr Dimensions cubed() const { return power(3); }

  constexpr Dimensions fractional_power(value_type d) const {
    return UNITS_ASSERT(exponents_divisible_by(d),
                        DimensionError("Cannot take root-" + std::to_string(d) +
                                       " of Dimensions " + to_string(true))),
           Dimensions(length() / d, mass() / d, time() / d, current() / d,
                      temperature() / d, luminosity() / d, amount() / d);
  }

  constexpr Dimensions sqrt() const { return fractional_power(2); }

  constexpr Dimensions cbrt() const { return fractional_power(3); }

  constexpr const value_type& length() const { return exponents_[kLength]; }
  constexpr const value_type& mass() const { return exponents_[kMass]; }
  constexpr const value_type& time() const { return exponents_[kTime]; }
  constexpr const value_type& current() const { return exponents_[kCurrent]; }
  constexpr const value_type& temperature() const {
    return exponents_[kTemperature];
  }
  constexpr const value_type& luminosity() const {
    return exponents_[kLuminosity];
  }
  constexpr const value_type& amount() const { return exponents_[kAmount]; }

  std::string to_string(bool write_if_none = false) const {
    return to_string({"l", "m", "t", "I", "T", "Iv", "N"}, write_if_none);
  }

  typedef std::array<const char*, kBaseDimensionCount> SymbolArray;

  std::string to_string(const SymbolArray& symbols,
                        bool write_if_none = false) const {
    std::stringstream stream;
    auto write_unit = [&stream](const std::string& symbol, int exponent,
                                bool with_space = true) {
      if (exponent) {
        if (with_space) {
          stream << " ";
        }
        stream << symbol;
        if (exponent != 1) {
          stream << exponent;
        }
        return true;
      }
      return false;
    };
    bool wrote = false;
    wrote |= write_unit(symbols[kLength], exponents_[kLength], wrote);
    wrote |= write_unit(symbols[kMass], exponents_[kMass], wrote);
    wrote |= write_unit(symbols[kTime], exponents_[kTime], wrote);
    wrote |= write_unit(symbols[kCurrent], exponents_[kCurrent], wrote);
    wrote |= write_unit(symbols[kTemperature], exponents_[kTemperature], wrote);
    wrote |= write_unit(symbols[kLuminosity], exponents_[kLuminosity], wrote);
    wrote |= write_unit(symbols[kAmount], exponents_[kAmount], wrote);
    if (!wrote && write_if_none) {
      stream << "<none>";
    }
    return stream.str();
  }

 private:
  constexpr tuple_type as_tuple() const {
    return std::make_tuple(length(), mass(), time(), current(), temperature(),
                           luminosity(), amount());
  }

  constexpr bool exponents_divisible_by(value_type n) const {
    return length() % n == 0 && mass() % n == 0 && time() % n == 0 &&
           current() % n == 0 && temperature() % n == 0 &&
           luminosity() % n == 0 && amount() % n == 0;
  }

  explicit constexpr Dimensions(const value_type& length_,
                                const value_type& mass_ = 0,
                                const value_type& time_ = 0,
                                const value_type& current_ = 0,
                                const value_type& temperature_ = 0,
                                const value_type& luminosity_ = 0,
                                const value_type& amount_ = 0)
      : exponents_({length_, mass_, time_, current_, temperature_, luminosity_,
                    amount_}) {}

  array_type exponents_;
  value_type padding_byte_ = 0;
};

inline Dimensions operator*(const Dimensions& lhs, const Dimensions& rhs) {
  Dimensions tmp(lhs);
  return tmp *= rhs;
}

inline Dimensions operator/(const Dimensions& lhs, const Dimensions& rhs) {
  Dimensions tmp(lhs);
  return tmp /= rhs;
}

inline Dimensions pow(const Dimensions& dims, Dimensions::value_type exponent) {
  return dims.power(exponent);
}

inline Dimensions sqrt(const Dimensions& dims) { return dims.sqrt(); }

inline Dimensions cbrt(const Dimensions& dims) { return dims.cbrt(); }

inline std::ostream& operator<<(std::ostream& stream, const Dimensions& dims) {
  stream << dims.to_string();
  return stream;
}

static_assert(std::is_standard_layout<Dimensions>::value,
              "Internal error: Dimensions is not standard layout");

template <typename T>
class Quantity {
 public:
  typedef T value_type;

  static value_type pi() { return T(3.1415926535897932384626L); }
  static value_type tau() { return T(2) * pi(); }

  // SI base units.
  static Quantity number() { return Quantity(Dimensions::None()); }
  static Quantity meter() { return Quantity(Dimensions::Length()); }
  static Quantity kilogram() { return Quantity(Dimensions::Mass()); }
  static Quantity second() { return Quantity(Dimensions::Time()); }
  static Quantity ampere() { return Quantity(Dimensions::Current()); }
  static Quantity kelvin() { return Quantity(Dimensions::Temperature()); }
  static Quantity candela() { return Quantity(Dimensions::Luminosity()); }
  static Quantity mole() { return Quantity(Dimensions::Amount()); }

  // SI derived units.
  static Quantity hertz() { return second().reciprocal(); }
  static Quantity radian() { return meter() / meter(); }
  static Quantity steradian() { return radian().squared(); }
  static Quantity newton() { return kilogram() * meter() / second().squared(); }
  static Quantity pascal() { return newton() / meter().squared(); }
  static Quantity joule() { return newton() * meter(); }
  static Quantity watt() { return joule() / second(); }
  static Quantity coulomb() { return ampere() * second(); }
  static Quantity volt() { return watt() / second(); }
  static Quantity farad() { return coulomb() / volt(); }
  static Quantity ohm() { return volt() / ampere(); }
  static Quantity siemens() { return ohm().reciprocal(); }
  static Quantity weber() { return joule() / ampere(); }
  static Quantity tesla() { return weber() / meter().squared(); }
  static Quantity henry() { return ohm() * second(); }
  static Quantity degree_celsius() { return kelvin(); }
  static Quantity lumen() { return candela() / steradian(); }
  static Quantity lux() { return lumen() / meter().squared(); }
  static Quantity becquerel() { return second().reciprocal(); }
  static Quantity gray() { return joule() / kilogram(); }
  static Quantity sievert() { return joule() / kilogram(); }
  static Quantity katal() { return mole() / second(); }

  // SI accepted units.
  static Quantity minute() { return T(60) * second(); }
  static Quantity hour() { return T(60) * minute(); }
  static Quantity day() { return T(24) * hour(); }
  static Quantity astronomical_unit() { return T(149597870700) * meter(); }
  static Quantity degree() { return pi() / T(180) * radian(); }
  static Quantity arc_minute() { return degree() / T(60); }
  static Quantity arc_second() { return arc_minute() / T(60); }
  static Quantity hectare() { return T(1e4) * meter().squared(); }
  static Quantity liter() { return T(1e-3) * meter().cubed(); }
  static Quantity metric_tonne() { return T(1e3) * kilogram(); }
  static Quantity dalton() { return T(1.660539040e-27) * kilogram(); }
  static Quantity electronvolt() { return T(1.602176634e-19) * joule(); }
  static Quantity speed_of_light() { return T(299792458) * meter() / second(); }
  static Quantity reduced_planck_constant() {
    return T(1.05457168e-34) * joule() * second();
  }
  static Quantity electron_mass() { return T(9.1093826e-31) * kilogram(); }
  static Quantity elementary_charge() { return T(1.60217653e-19) * coulomb(); }
  static Quantity planck_time() { return T(1.2880886677e-21) * second(); }
  static Quantity bohr_radius() { return T(5.291772108e-11) * meter(); }
  static Quantity hartree_energy() { return T(4.35974417e-18) * joule(); }

  // Unofficial units.
  static Quantity angstrom() { return T(1e-10) * meter(); }
  static Quantity are() { return T(100) * meter().squared(); }
  static Quantity barn() { return T(1e-28) * meter().squared(); }
  static Quantity bar() { return T(1e5) * pascal(); }
  static Quantity atmosphere() { return T(101325) * pascal(); }
  static Quantity barye() { return T(0.1) * pascal(); }
  static Quantity millimeter_of_mercury() {
    return T(133.322387415) * pascal();
  }
  static Quantity conventional_mercury() {
    return T(13595.1) * standard_gravity() * kilogram() / meter().cubed();
  }
  static Quantity torr() { return atmosphere() / T(760); }
  static Quantity year() { return T(365.25) * day(); }

  // CGS units.
  static Quantity gal() { return T(1e-2) * meter() / second().squared(); }
  static Quantity dyne() { return T(1e-5) * newton(); }
  static Quantity erg() { return T(1e-7) * joule(); }
  static Quantity poise() { return T(0.1) * pascal() * second(); }
  static Quantity stokes() { return T(1e-4) * meter().squared() / second(); }
  static Quantity gauss() { return T(1e-4) * tesla(); }
  static Quantity statcoulomb() { return T(3.33564e-10) * coulomb(); }
  static Quantity franklin() { return statcoulomb(); }
  static Quantity statvolt() { return T(299.792458) * volt(); }
  static Quantity oersted() { return T(79.57747) * ampere() / meter(); }
  static Quantity phot() { return T(1e4) * lux(); }
  static Quantity rad() { return T(1e-2) * gray(); }
  static Quantity roentgen_equivalent_man() { return T(1e-2) * sievert(); }
  static Quantity maxwell() { return T(1e-8) * weber(); }

  // Misc units.
  static Quantity week() { return T(7) * day(); }
  static Quantity pounds_per_square_inch() { return T(6.894757e3) * pascal(); }
  static Quantity nutritional_calorie() { return T(4.186e3) * joule(); }
  static Quantity curie() { return T(3.7e10) * becquerel(); }
  static Quantity dobson_unit() {
    return T(0.4462e-3) * mole() / meter().squared();
  }
  static Quantity furlong() { return T(201.168) * meter(); }
  static Quantity horsepower() { return T(745.7) * watt(); }
  static Quantity roentgen() {
    return T(2.58) * ampere() * second() / kilogram();
  }
  static Quantity rack_unit() { return T(44.45e-3) * meter(); }
  static Quantity turn() { return T(2) * pi() * radian(); }
  static Quantity british_thermal_units() {
    return T(1.05505585262e3) * joule();
  }
  static Quantity revolutions_per_minute() { return turn() / minute(); }
  static Quantity counts_per_second() { return second().reciprocal(); }
  static Quantity pound() { return T(0.45359237) * kilogram(); }
  static Quantity pound_of_force() { return T(4.448222) * newton(); }
  static Quantity tex() { return T(1e-6) * kilogram() / meter(); }
  static Quantity denier() { return tex() / T(9); }
  static Quantity rutherford() { return T(1e6) * becquerel(); }
  static Quantity eon() { return T(1e9) * year(); }
  static Quantity inch() { return T(25.4e-3) * meter(); }
  static Quantity foot() { return T(12) * inch(); }
  static Quantity yard() { return T(3) * foot(); }
  static Quantity mile() { return T(1760) * yard(); }
  static Quantity standard_gravity() {
    return T(9.806650) * meter() / second().squared();
  }
  static Quantity vacuum_permeability() {
    return T(4e-7) * pi() * newton() / ampere().squared();
  }
  static Quantity vacuum_permittivity() {
    return (vacuum_permeability() * speed_of_light().squared()).reciprocal();
  }
  static Quantity vacuum_impedance() {
    return vacuum_permeability() * speed_of_light();
  }
  static Quantity gravitational_constant() {
    return T(6.67408e-11) * meter().cubed() / kilogram() / second().squared();
  }
  static Quantity avogadro_number() { return T(6.02214129e23) / mole(); }
  static Quantity boltzmann_constant() {
    return T(1.3806488e-23) * joule() / kelvin();
  }
  static Quantity stefan_boltzmann_constant() {
    return T(5.670373e-8) * watt() / meter().squared() /
           kelvin().squared().squared();
  }
  static Quantity neutron_mass() { return T(1.674927351e-27) / kilogram(); }
  static Quantity proton_mass() { return T(1.672621777e-27) / kilogram(); }

  // Astronomy units.
  static Quantity solar_mass() { return T(1.98847e30) * kilogram(); }
  static Quantity earth_mass() { return T(5.9722e24) * kilogram(); }
  static Quantity jupiter_mass() { return T(1.898130e27) * kilogram(); }
  static Quantity solar_radius() { return T(6.957e8) * meter(); }
  static Quantity earth_radius() { return T(6.3781e6) * meter(); }
  static Quantity jupiter_radius() { return T(7.1492e7) * meter(); }
  static Quantity light_year() { return T(9460730472580800) * meter(); }
  static Quantity parsec() { return T(648000) / pi() * astronomical_unit(); }
  static Quantity lunar_distance() { return T(3.84402e8) * meter(); }
  static Quantity solar_luminosity() { return T(3.939e26) * watt(); }
  static Quantity jansky() {
    return T(1e-26) * watt() / meter().squared() / hertz();
  }

  typedef std::unordered_map<std::string, Quantity> symbol_map_type;

  // TODO: Consider adding long-form names.
  //{"meter", meter()},
  //{"metre", meter()},
  // These form a self-consistent grammar where any unit can be combined with
  // any prefix without colliding with another unit name.
  static const symbol_map_type& unit_symbol_map() {
    static const symbol_map_type units = {
        // SI base units.
        {"m", meter()},
        {"g", T(1e-3) * kilogram()},
        {"s", second()},
        {"A", ampere()},
        {"K", kelvin()},
        {"cd", candela()},
        {"mol", mole()},
        // SI derived units.
        {"Hz", hertz()},
        {"rad", radian()},
        {"sr", steradian()},
        {"N", newton()},
        {"Pa", pascal()},
        {"J", joule()},
        {"W", watt()},
        {"C", coulomb()},
        {"V", volt()},
        {"F", farad()},
        {"ohm", ohm()},
        {"S", siemens()},
        {"Wb", weber()},
        {"T", tesla()},
        {"H", henry()},
        {"degC", degree_celsius()},
        {"\xB0"  // degree symbol
         "C",
         degree_celsius()},
        {"lm", lumen()},
        {"lx", lux()},
        {"Bq", becquerel()},
        {"Gy", gray()},
        {"Sv", sievert()},
        {"kat", katal()},
        // SI acceptable units.
        {"min", minute()},
        {"h", hour()},
        {"day", day()},  // Can't use "d" because of candela vs. centi-days.
        {"deg", degree()},
        {"\xB0", degree()},  // degree symbol
        {"'", arc_minute()},
        {"\"", arc_second()},
        {"ha", hectare()},
        {"L", liter()},
        {"l", liter()},
        {"t", metric_tonne()},
        {"AU", astronomical_unit()},
        {"eV", electronvolt()},
        {"Da", dalton()},
        {"u", dalton()},
        {"c_0", speed_of_light()},
        {"h_", reduced_planck_constant()},
        {"m_e", electron_mass()},
        {"e", elementary_charge()},
        {"a_0", bohr_radius()},
        {"E_h", hartree_energy()},
        // Unofficial units.
        {"Angstrom", angstrom()},
        {"\u212B", angstrom()},  // angstrom symbol
        {"\xC5", angstrom()},    // letter A with circle over it
        {"are", are()},
        {"b", barn()},
        {"bar", bar()},
        {"atm", atmosphere()},
        {"Ba", barye()},
        {"Hg", conventional_mercury()},
        {"mmHg", millimeter_of_mercury()},
        {"Torr", torr()},
        {"year", year()},
        {"yr", year()},
        // Astronomy units.
        {"Msun", solar_mass()},
        {"M_S", solar_mass()},
        {"M\u2609", solar_mass()},  // circle with dot in it
        {"Mearth", earth_mass()},
        {"M_E", earth_mass()},
        {"M\u2295", earth_mass()},  // circle with plus in it
        {"Mjupiter", jupiter_mass()},
        {"M_J", jupiter_mass()},
        {"Rsun", solar_radius()},
        {"R_S", solar_radius()},
        {"R\u2609", solar_radius()},  // circle with dot in it
        {"Rearth", earth_radius()},
        {"R_E", earth_radius()},
        {"R\u2295", earth_radius()},  // circle with plus in it
        {"Rjupiter", jupiter_radius()},
        {"R_J", jupiter_radius()},
        {"ly", light_year()},
        {"pc", parsec()},
        {"LD", lunar_distance()},
        {"Lsun", solar_luminosity()},
        {"Jy", jansky()},
        // CGS units.
        {"Gal", gal()},
        {"dyn", dyne()},
        {"erg", erg()},
        {"P", poise()},
        {"St", stokes()},
        {"gauss", gauss()},
        {"Fr", franklin()},
        {"esu", franklin()},
        {"statV", statvolt()},
        {"Oe", oersted()},
        {"phot", phot()},
        {"rd", rad()},
        {"rem", roentgen_equivalent_man()},
        {"Mx", maxwell()},
        // Misc units.
        {"wk", week()},
        {"psi", pounds_per_square_inch()},
        {"Cal", nutritional_calorie()},
        {"Ci", curie()},
        {"DU", dobson_unit()},
        {"fur", furlong()},
        {"hp", horsepower()},
        {"R", roentgen()},
        {"U", rack_unit()},
        {"tr", turn()},
        {"rev", turn()},
        {"cyc", turn()},
        {"BTU", british_thermal_units()},
        {"Btu", british_thermal_units()},
        {"rpm", revolutions_per_minute()},
        {"cps", counts_per_second()},
        {"lb", pound()},
        {"lbf", pound_of_force()},
        {"tex", tex()},
        {"D", denier()},
        {"Rd", rutherford()},
        {"eon", eon()},
        {"inch", inch()},
        {"foot", foot()},
        {"yard", yard()},
        {"yd", yard()},
        {"mile", mile()},
        {"mi", mile()},
        {"g_0", standard_gravity()},
        {"mu_0", vacuum_permeability()},
        {"epsilon_0", vacuum_permittivity()},
        {"Z_0", vacuum_impedance()},
        {"N_A", avogadro_number()},
        {"k", boltzmann_constant()},
        {"\u03C3", stefan_boltzmann_constant()},  // lowercase sigma
        {"m_n", neutron_mass()},
        {"m_p", proton_mass()},
        {"G", gravitational_constant()},
    };
    return units;
  }

  // SI prefixes.
  static value_type yotta() { return T(1e24); }
  static value_type zetta() { return T(1e21); }
  static value_type exa() { return T(1e18); }
  static value_type peta() { return T(1e15); }
  static value_type tera() { return T(1e12); }
  static value_type giga() { return T(1e9); }
  static value_type mega() { return T(1e6); }
  static value_type kilo() { return T(1e3); }
  static value_type hecto() { return T(1e2); }
  static value_type deka() { return T(1e1); }
  static value_type deci() { return T(1e-1); }
  static value_type centi() { return T(1e-2); }
  static value_type milli() { return T(1e-3); }
  static value_type micro() { return T(1e-6); }
  static value_type nano() { return T(1e-9); }
  static value_type pico() { return T(1e-12); }
  static value_type femto() { return T(1e-15); }
  static value_type atto() { return T(1e-18); }
  static value_type zepto() { return T(1e-21); }
  static value_type yocto() { return T(1e-24); }

  // Binary prefixes.
  static value_type kibi() { return T(1024); }
  static value_type mebi() { return T(1024 * kibi()); }
  static value_type gibi() { return T(1024 * mebi()); }
  static value_type tebi() { return T(1024 * gibi()); }
  static value_type pebi() { return T(1024 * tebi()); }
  static value_type exbi() { return T(1024 * pebi()); }
  static value_type zebi() { return T(1024 * exbi()); }
  static value_type yobi() { return T(1024 * zebi()); }

  // TODO: Consider making these value_type instead of Quantity.
  static const symbol_map_type& prefix_symbol_map() {
    static const symbol_map_type prefixes = {
        {"Y", yotta()}, {"Z", zetta()}, {"E", exa()},      {"P", peta()},
        {"T", tera()},  {"G", giga()},  {"M", mega()},     {"k", kilo()},
        {"h", hecto()}, {"da", deka()}, {"d", deci()},     {"c", centi()},
        {"m", milli()}, {"u", micro()}, {"\xB5", micro()}, {"n", nano()},
        {"p", pico()},  {"f", femto()}, {"a", atto()},     {"z", zepto()},
        {"y", yocto()}, {"Ki", kibi()}, {"Mi", mebi()},    {"Gi", gibi()},
        {"Ti", tebi()}, {"Pi", pebi()}, {"Ei", exbi()},    {"Zi", zebi()},
        {"Yi", yobi()},
    };
    return prefixes;
  }

  Quantity() : value_() {}

  // Allow implicit conversion from value_type (dimensionless).
  Quantity(const value_type& value) : value_(value), dims_() {}

  // Allow implicit conversion from Dimensions (unit magnitude).
  Quantity(const Dimensions& dims) : value_(1), dims_(dims) {}

  Quantity(const value_type& value, const Dimensions& dims)
      : value_(value), dims_(dims) {}

  explicit Quantity(const std::string& units) {
    bool units_parsed_successfully = parse_units(units, this);
    assert(units_parsed_successfully);
    (void)units_parsed_successfully;
  }

  Quantity(const value_type& value, const std::string& units)
      : Quantity(value) {
    Quantity q;
    bool units_parsed_successfully = parse_units(units, &q);
    assert(units_parsed_successfully);
    (void)units_parsed_successfully;
    *this *= q;
  }

  const value_type& magnitude() const { return value_; }

  value_type magnitude(const Quantity& in_units) const {
    UNITS_ASSERT(
        in_units.dimensions() == dimensions(),
        QuantityError("Dimension mismatch: " + dimensions().to_string(true) +
                      " vs. " + in_units.dimensions().to_string(true)));
    return static_cast<value_type>(*this / in_units);
  }

  value_type magnitude(const std::string& in_units) const {
    return magnitude(Quantity(in_units));
  }

  const Dimensions& dimensions() const { return dims_; }

  bool operator==(const Quantity& rhs) const {
    return dims_ == rhs.dims_ && value_ == rhs.value_;
  }

  bool operator!=(const Quantity& rhs) const {
    return dims_ != rhs.dims_ || value_ != rhs.value_;
  }

  bool operator<(const Quantity& rhs) const {
    // TODO: Factor out this assertion into a helper method.
    UNITS_ASSERT(
        rhs.dimensions() == dimensions(),
        QuantityError("Dimension mismatch: " + dimensions().to_string(true) +
                      " vs. " + rhs.dimensions().to_string(true)));
    return value_ < rhs.value_;
  }

  bool operator<=(const Quantity& rhs) const {
    return *this == rhs || *this < rhs;
  }

  bool operator>(const Quantity& rhs) const { return !(*this <= rhs); }

  bool operator>=(const Quantity& rhs) const { return !(*this < rhs); }

  Quantity operator-() const { return Quantity(-value_, dims_); }

  Quantity& operator*=(const Quantity& rhs) {
    dims_ *= rhs.dims_;
    value_ *= rhs.value_;
    return *this;
  }

  Quantity& operator/=(const Quantity& rhs) {
    dims_ /= rhs.dims_;
    value_ /= rhs.value_;
    return *this;
  }

  Quantity& operator+=(const Quantity& rhs) {
    UNITS_ASSERT(dims_ == rhs.dims_,
                 QuantityError("Dimension mismatch: " + dims_.to_string(true) +
                               " vs. " + rhs.dims_.to_string(true)));
    value_ += rhs.value_;
    return *this;
  }

  Quantity& operator-=(const Quantity& rhs) {
    UNITS_ASSERT(dims_ == rhs.dims_,
                 QuantityError("Dimension mismatch: " + dims_.to_string(true) +
                               " vs. " + rhs.dims_.to_string(true)));
    value_ -= rhs.value_;
    return *this;
  }

  // Dimensionless value.
  explicit operator value_type() const {
    UNITS_ASSERT(!dims_, QuantityError("Quantity is not dimensionless: " +
                                       dims_.to_string(true)));
    return value_;
  }

  Quantity power(value_type n) const {
    // Only integer exponents are supported.
    Dimensions::value_type ni(n);
    UNITS_ASSERT(ni == n, QuantityError("Exponent is not representable: " +
                                        std::to_string(n)));
    return Quantity(std::pow(value_, n), dims_.power(n));
  }

  Quantity fractional_power(Dimensions::value_type d) const {
    return Quantity(std::pow(value_, value_type(1) / d),
                    dims_.fractional_power(d));
  }

  Quantity reciprocal() const { return power(-1); }

  Quantity squared() const { return power(2); }

  Quantity cubed() const { return power(3); }

  Quantity sqrt() const { return Quantity(::sqrt(value_), dims_.sqrt()); }

  Quantity cbrt() const { return Quantity(::cbrt(value_), dims_.cbrt()); }

  Quantity units() const { return Quantity(dimensions()); }

 private:
  value_type value_;
  Dimensions dims_;
};

static_assert(std::is_standard_layout<Quantity<double>>::value,
              "Internal error: Quantity is not standard layout");

#define DEFINE_QUANTITY_BINARY_OPERATOR(op)                              \
  template <typename T>                                                  \
  inline Quantity<T> operator op(const Quantity<T>& lhs,                 \
                                 const Quantity<T>& rhs) {               \
    Quantity<T> tmp(lhs);                                                \
    return tmp op## = rhs;                                               \
  }                                                                      \
  template <typename T>                                                  \
  inline Quantity<T> operator op(const T& lhs, const Quantity<T>& rhs) { \
    Quantity<T> tmp(lhs);                                                \
    return tmp op## = rhs;                                               \
  }                                                                      \
  template <typename T>                                                  \
  inline Quantity<T> operator op(const Quantity<T>& lhs, const T& rhs) { \
    Quantity<T> tmp(lhs);                                                \
    return tmp op## = Quantity<T>(rhs);                                  \
  }

DEFINE_QUANTITY_BINARY_OPERATOR(*)
DEFINE_QUANTITY_BINARY_OPERATOR(/)
DEFINE_QUANTITY_BINARY_OPERATOR(+)
DEFINE_QUANTITY_BINARY_OPERATOR(-)

#undef DEFINE_QUANTITY_BINARY_OPERATOR

template <typename T>
inline Quantity<T> pow(const Quantity<T>& quantity,
                       typename Quantity<T>::value_type exponent) {
  return quantity.power(exponent);
}

template <typename T>
inline Quantity<T> fabs(const Quantity<T>& quantity) {
  return Quantity<T>(::fabs(quantity.magnitude()), quantity.dimensions());
}

template <typename T>
inline Quantity<T> sqrt(const Quantity<T>& quantity) {
  return quantity.sqrt();
}

template <typename T>
inline Quantity<T> cbrt(const Quantity<T>& quantity) {
  return quantity.cbrt();
}

template <typename T>
inline std::ostream& operator<<(std::ostream& stream, const Quantity<T>& q) {
  stream << q.magnitude();
  if (q.dimensions()) {
    stream << " "
           << q.dimensions().to_string({"m", "kg", "s", "A", "K", "cd", "mol"});
  }
  return stream;
}

using Quantity64 = Quantity<double>;

namespace detail {

template <typename T>
inline std::string make_reverse_regex_from_keys(
    const typename Quantity<T>::symbol_map_type& symbol_map) {
  // Keys are reverse-sorted so that greedy regex parsing matches long symbols
  // first before short ones (which may be prefixes of long ones).
  std::vector<std::string> rsorted_symbols;
  rsorted_symbols.reserve(symbol_map.size());
  for (const auto& symbol_quantity : symbol_map) {
    const std::string& symbol = symbol_quantity.first;
    std::string rsymbol(symbol.rbegin(), symbol.rend());
    rsorted_symbols.push_back(rsymbol);
  }
  std::sort(rsorted_symbols.rbegin(), rsorted_symbols.rend());

  std::stringstream ss;
  ss << "(";
  bool is_first_symbol = true;
  for (const std::string& symbol : rsorted_symbols) {
    if (is_first_symbol) {
      is_first_symbol = false;
    } else {
      ss << "|";
    }
    ss << symbol;
  }
  ss << ")";
  return ss.str();
}

inline std::string reverse(const std::string& str) {
  return std::string(str.rbegin(), str.rend());
}

inline bool is_whitespace(const std::string& str) {
  return std::find_if(str.begin(), str.end(),
                      [](char c) { return !std::isspace(c); }) == str.end();
}

inline bool contains_nonseparator_char(const std::string& str) {
  return std::find_if(str.begin(), str.end(), [](char c) {
           const char middot = '\xB7';
           const char times = '\xD7';
           return !(std::isspace(c) || c == '*' || c == ',' || c == middot ||
                    c == times);
         }) != str.end();
}

}  // end namespace detail

// Returns false on error, or throws an exception if RTUNITS_USE_EXCEPTIONS=1.
// Valid separator symbols are: whitespace, asterisk, comma, middot, times.
template <typename T>
inline bool parse_units(const std::string& str, Quantity<T>* result_ptr) {
  // TODO: Add support for a numerical scale factor prefix.
  // This uses a regex to match each unit with the form:
  //   [division_symbol] [prefix] unit [exponent_symbol] [exponent]
  // Whitespace within this form is optional, but units must be separated by
  // a separator symbol (or a division symbol).
  // Note that the regex is applied in reverse so that units are matched first
  // before prefixes.
  static const std::regex re = []() {
    const std::string div_symbol = R"((\/)?)";
    const std::string ws = R"(\s*)";  // whitespace
    const std::string exp_symbol = R"((\*\*|\^)?)";
    const std::string exponent = R"((\d+-?)*)";
    const std::string prefix = detail::make_reverse_regex_from_keys<T>(
                                   Quantity<T>::prefix_symbol_map()) +
                               "?";
    // Positive lookahead to ensure that the next symbol is a separator (or
    // division symbol).
    // \xB7 = middle dot, \xD7 == multiplication sign.
    const std::string unit_boundary = R"((?=[\s*,\/\xB7\xD7]|$))";
    const std::string unit =
        detail::make_reverse_regex_from_keys<T>(Quantity<T>::unit_symbol_map());
    return std::regex(exponent + ws + exp_symbol + ws + unit + prefix +
                      unit_boundary + ws + div_symbol);
  }();
  if (detail::is_whitespace(str)) {
    *result_ptr = Quantity<T>(Dimensions::None());
    return true;
  }
  using regex_iter = std::sregex_iterator;
  Quantity<T> result(1);
  std::string rstr = detail::reverse(str);
  auto it = regex_iter(rstr.begin(), rstr.end(), re);
  bool any_units_matched = it != regex_iter();
  if (!any_units_matched) {
    return UNITS_THROW(QuantityParseError("No units matched in string: \"" +
                                          str + "\"")),
           false;
  }
  regex_iter last_it;
  for (; it != regex_iter(); ++it) {
    last_it = it;
    const std::smatch& match = *it;
    if (detail::contains_nonseparator_char(match.prefix().str())) {
      std::string fail_str = detail::reverse(match.prefix().str());
      return UNITS_THROW(
                 QuantityParseError("Could not parse \"" + fail_str + "\"")),
             false;
    }
    // Note that the first submatch is the whole match.
    bool valid_unit_match = match.size() == 6;
    if (!valid_unit_match) {
      return UNITS_THROW(QuantityParseError("Bad unit match: " + match.str())),
             false;
    }
    int m0 = 1;
    std::string div_symbol = detail::reverse(match[m0 + 4]);
    std::string prefix = detail::reverse(match[m0 + 3]);
    std::string unit = detail::reverse(match[m0 + 2]);
    std::string exponent_str = detail::reverse(match[m0 + 0]);
    Quantity<T> term(1);
    term *= Quantity<T>::unit_symbol_map().at(unit);
    if (!prefix.empty()) {
      term *= Quantity<T>::prefix_symbol_map().at(prefix);
    }
    if (!exponent_str.empty()) {
      term = term.power(std::atoi(exponent_str.c_str()));
    }
    if (!div_symbol.empty()) {
      term = term.reciprocal();
    }
    result *= term;
  }
  if (detail::contains_nonseparator_char(last_it->suffix().str())) {
    std::string fail_str = detail::reverse(last_it->suffix().str());
    return UNITS_THROW(
               QuantityParseError("Could not parse \"" + fail_str + "\"")),
           false;
  }
  *result_ptr = result;
  return true;
}

// Returns false on error, or throws an exception if RTUNITS_USE_EXCEPTIONS=1.
template <typename T>
inline bool convert_units(T* magnitude, const std::string& from_units,
                          const std::string& to_units) {
  Quantity<T> from_q;
  if (!parse_units(from_units, &from_q)) {
    return false;
  }
  Quantity<T> to_q;
  if (!parse_units(to_units, &to_q)) {
    return false;
  }
  if (from_q.dimensions() != to_q.dimensions()) {
    return UNITS_THROW(QuantityError(
               "Dimension mismatch: " + from_q.dimensions().to_string(true) +
               " vs. " + to_q.dimensions().to_string(true))),
           false;
  }
  *magnitude *= static_cast<T>(from_q / to_q);
  return true;
}

namespace detail {

// As used in Boost.
template <class T>
inline void hash_combine(std::size_t& seed, const T& v) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

}  // end namespace detail

}  // end namespace rtunits

namespace std {

template <>
struct hash<::rtunits::Dimensions> {
  size_t operator()(const ::rtunits::Dimensions& dims) const {
    using value_type = ::rtunits::Dimensions::value_type;
    size_t seed = 0;
    for (const auto& e : dims.exponents()) {
      ::rtunits::detail::hash_combine(seed, hash<value_type>()(e));
    }
    return seed;
  }
};

template <typename T>
struct hash<rtunits::Quantity<T>> {
  size_t operator()(const ::rtunits::Quantity<T>& q) const {
    size_t seed = hash<T>()(q.magnitude());
    ::rtunits::detail::hash_combine(
        seed, hash<::rtunits::Dimensions>()(q.dimensions()));
    return seed;
  }
};

}  // end namespace std

#undef UNITS_THROW
#undef UNITS_ASSERT
#undef CONSTEXPR_ASSERT
