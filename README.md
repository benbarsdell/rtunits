# RT-Units

[![Build Status](https://travis-ci.com/benbarsdell/rtunits.svg?branch=master)](https://travis-ci.com/benbarsdell/rtunits)

RT-Units - A single-header C++11 library for working with physical
quantities at runtime, including support for parsing and
serialization.

The `rtunits::Quantity` template represents a physical quantity with a
magnitude and dimensions. It supports an extensive set of physical
units and scale prefixes, and implements robust parsing and
serialization.

# Example

```C++
#include <rtunits.hpp>
#include <iostream>

using rtunits::Quantity64;

int main(int argc, char* argv[]) {
  Quantity64 force(1.25, "kN");
  Quantity64 displacement(12.0, "mm");
  Quantity64 work = force * displacement + 5. * Quantity64::joule();
  std::cout << work << std::endl;
  std::cout << (work == Quantity64(20, "J")) << std::endl;
  std::cout << (work.magnitude("mW s") == 20e3) << std::endl;
  return 0;
}
```

Output:

```
20 m2 kg s-2
1
1
```

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
