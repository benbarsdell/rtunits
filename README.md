# RT-Units

[![Build Status](https://github.com/benbarsdell/rtunits/actions/workflows/ci.yml/badge.svg?event=push)](https://github.com/benbarsdell/rtunits/actions?query=workflow%3ACI)


RT-Units - A single-header C++11 library for working with physical
quantities at runtime, including support for parsing and
serialization.

The `rtunits::Quantity` template represents a physical quantity with a
magnitude and dimensions. It supports an extensive set of physical
units and scale prefixes, and implements robust parsing and
serialization.

# Example

```C++
#include <iostream>
#include <rtunits.hpp>

using rtunits::Quantity64;
using rtunits::SpecificQuantity64;
using Dims = rtunits::Dimensions;

int main(int argc, char* argv[]) {
  Quantity64 force(1.25, "kN");
  Quantity64 displacement(12.0, "mm");
  Quantity64 work = force * displacement + 5. * Quantity64::joule();
  std::cout << "Work: " << work << '\n';
  std::cout << "Dims: " << work.dimensions() << '\n';
  std::cout << (work.dimensions() == Dims::Force() * Dims::Length()) << '\n';
  std::cout << (work == Quantity64(20, "J")) << '\n';
  std::cout << (work.value("mW s") == 20e3) << '\n';

  SpecificQuantity64 dm(500., "pc/cm^3");
  dm *= 2;
  std::cout << "DM: " << dm << '\n';
  std::cout << "  = " << static_cast<Quantity64>(dm) << std::endl;
}
```

Output:

```
Work: 20 m^2 kg s^-2
Dims: L^2 M T^-2
1
1
1
DM: 1000 pc/cm^3
  = 3.08568e+25 m^-2
```

# Performance

Units are parsed using a very large regular expression, which can be slow.
To improve performance, consider doing the following:

- Define RTUNITS_USE_BOOST_REGEX=1 to use [boost::regex](https://github.com/boostorg/regex),
  which is 4-8x faster than std::regex. (Note that boost::regex can be used as a
  standalone header-only library without requiring the rest of Boost).
- Define RTUNITS_NO_ASTRONOMY_UNITS=1, RTUNITS_NO_CGS_UNITS=1, and/or
  RTUNITS_NO_MISC_UNITS=1 before including `rtunits.hpp` to disable support for
  unneeded units. This can speed up parsing by up to 2x.


# Error handling

Exceptions are used by default, but can be disabled by defining
`RTUNITS_USE_EXCEPTIONS` to 0 before including the rtunits header:

```
#define RTUNITS_USE_EXCEPTIONS 0
#include <rtunits.hpp>
```

# Tests

Run tests using the following commands:

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make check
