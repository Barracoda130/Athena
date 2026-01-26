#include "Athena.hpp"
#include "MpfrPlayground.hpp"
#include <iostream>

int main()
{
	Athena::Number num1( "123", 256 );
	Athena::Number num2( "9999999999999999999", 256 );

	std::cout << "num1: " << num1 << std::endl;
	std::cout << "num2: " << num2 << std::endl;

	return 0;
}