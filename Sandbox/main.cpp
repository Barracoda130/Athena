#include "Athena.hpp"
#include "MpfrPlayground.hpp"
#include <iostream>

int main()
{
	std::string a("123");
	std::string b( "456" );

	number_str_t n1 = toNumberStr( "123");
	number_str_t n2 = toNumberStr( "456");
	number_str_t n3 = strAdd( n1, n2 );

	std::cout << fromNumberStr( n3 ) << std::endl;

	return 0;
}