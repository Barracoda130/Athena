#include "Athena.hpp"
#include "MpfrPlayground.hpp"
#include <iostream>

int main()
{
	Athena::testFunction();

	//add();
	std::cout << sizeof( unsigned long long ) << std::endl;

	Athena::Number n1( 50, 100 );
	Athena::Number n2( 60, 100 );
	Athena::Number n3( 70, 100 );

	add( n3, n1, n2, Athena::RNDA );

	n1.set( n2 );

	std::cout << "end" << std::endl;

	return 0;
}