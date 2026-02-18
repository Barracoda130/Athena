#pragma once

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4146)
#endif

#include <mpfr.h>

#ifdef _MSC_VER
#pragma warning(pop)
#endif


std::ostream& operator<<( std::ostream& os, const mpfr_t a_Num );