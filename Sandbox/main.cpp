#include "Athena.hpp"
#include "MpfrPlayground.hpp"
#include <iostream>

int main()
{
	Athena::Number num1( "-12", 256 );
	Athena::Number num2( "13", 256 );

	Athena::Number result( 256 );
	add( result, num1, num2, MPFR_RNDN );
	std::cout << "result: " << result << std::endl;

	std::cout << "num1: " << num1 << std::endl;
	std::cout << "num2: " << num2 << std::endl;

	mpfr_t a, b, mpfr_r;
	mpfr_init2( a, 256 );
	mpfr_init2( b, 256 );
	mpfr_init2( mpfr_r, 256 );
	mpfr_set_str( a, "123", 10, MPFR_RNDN );
	mpfr_set_str( b, "9999999999999999999", 10, MPFR_RNDN );

	mpfr_add( mpfr_r, a, b, MPFR_RNDN );

	bool eq = static_cast<bool>( mpfr_equal_p( a, b ) );
	std::cout << "mpfr_equal_p result: " << eq << std::endl;
	std::cout << "eq " << ( mpfr_r == result ) << std::endl;
	std::cout << a << std::endl;

	add();

	// Default: 100 limbs, 10000 iterations
	//Athena::benchmark_add_n_implementations();

	// Custom parameters
	Athena::benchmark_add_n_implementations( 50, 100000 );
	return 0;

	return 0;
}