#include "MpfrInclude.hpp"

#include <iostream>
#include <format>

// Method to allow the easy printing of the mpfr type using C++ streams. This is used for testing and debugging purposes.
std::ostream& operator<<( std::ostream& os, mpfr_srcptr a_Num )
{
	if ( a_Num == nullptr )
		return os << "<null mpfr>";

	mpfr_exp_t exponent;
	char* mantissa = mpfr_get_str( nullptr, &exponent, 10, 0, a_Num, MPFR_RNDN );
	if ( mantissa == nullptr )
		return os << "<mpfr_get_str failed>";

	char initial = mantissa[0];
	char* digits = mantissa + 1;

	std::string out = std::format( "{}.{}e{}", initial, digits, exponent - 1 );
	mpfr_free_str( mantissa );

	return os << out;
}