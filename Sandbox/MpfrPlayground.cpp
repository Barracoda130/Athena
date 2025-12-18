#include "MpfrPlayground.hpp"


#include "MpfrInclude.hpp"
#include <cstdio>

void add()
{
	mpfr_t n1, n2;
	mpfr_init2( n1, 256 );
	mpfr_init2( n2, 256 );
	
	mpfr_set_str( n1, "1234567890987654321", 10, MPFR_RNDD );
	mpfr_set_str( n2, "9999999999999999999", 10, MPFR_RNDD );

	// Print using mpfr_printf (use %Rf for mpfr_t)
	mpfr_printf("n1 = %.0Rf\n", n1); // "%.0Rf" prints with zero fractional digits (integer style)
	mpfr_printf("n2 = %.0Rf\n", n2);

	// Do the addition (result stored in n1)
	mpfr_add( n1, n1, n2, MPFR_RNDU );

	// Print the result
	mpfr_printf("sum = %.0Rf\n", n1);

	// Alternative: print to a FILE* with mpfr_out_str (base 10, 0 = all digits)
	mpfr_out_str(stdout, 10, 0, n1, MPFR_RNDN);
	putchar('\n');

	mpfr_clear( n1 );
	mpfr_clear( n2 );
}