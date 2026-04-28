#pragma once

// Suppress warnings from mpfr.h, which is a C library and may have some constructs that trigger warnings in C++

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4146)
#endif

#include <mpfr.h>

#ifdef _MSC_VER
#pragma warning(pop)
#endif


std::ostream& operator<<( std::ostream& os, const mpfr_t a_Num );