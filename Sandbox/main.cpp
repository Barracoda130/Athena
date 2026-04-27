#include "Athena.hpp"
#include "MpfrPlayground.hpp"
#include <iostream>
#include <stdio.h>



void atnThings()
{
	Athena::Number n1( "12345678987654321", 100 );
	Athena::Number n2( "12345678987654321", 100 );
	Athena::Number result( 100 );
	add( result, n1, n2, MPFR_RNDN );

	Athena::Number n3( "12345678987654321", 200 );
	Athena::Number n4( "12345678987654321", 200 );
	Athena::Number result2( 200 );
	add( result2, n3, n4, MPFR_RNDN );

	Athena::Number n5( "12345678987654321", 300 );
	Athena::Number n6( "12345678987654321", 300 );
	Athena::Number result3( 300 );
	add( result3, n5, n6, MPFR_RNDN );

	Athena::Number n7( "12345678987654321", 400 );
	Athena::Number n8( "12345678987654321", 400 );
	Athena::Number result4( 400 );
	add( result4, n7, n8, MPFR_RNDN );
}

int main()
{
	atnThings();
	return 0;
}