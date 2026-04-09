#include "MpfrPlayground.hpp"


#include "MpfrInclude.hpp"
#include <cstdio>
#include <bitset>
#include <iostream>

void add()
{
	mpfr_t n1, n2;
	mpfr_init2( n1, 256 );
	mpfr_init2( n2, 256 );
	
	mpfr_set_str( n1, "123", 16, MPFR_RNDD );
	mpfr_set_str( n2, "ff1234567890192837ff", 16, MPFR_RNDD );
	unsigned long long b = 0x123456;
	unsigned long long* a = &b;
	

	std::bitset<256> bits_n1(n1->_mpfr_d);
	std::bitset<256> bits_n2(n2->_mpfr_d);	
	std::cout << "Bits of n1: " << bits_n1 << std::endl;
	std::cout << "Bits of n2: " << bits_n2 << std::endl;

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