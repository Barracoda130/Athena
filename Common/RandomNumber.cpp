#include "RandomNumber.hpp"

#include <random>

int randInt( int a_Low, int a_High )
{
	return (std::rand() % (a_High + 1 - a_Low)) + a_Low;
}

void generateInt( std::size_t a_Length, char* result )
{
	for ( index_t i = 0; i < a_Length; i++ )
	{
		result[i] = static_cast<char>( randInt( 0, 9 ) ) + '0';	// Get random number between 0 and 10
	}
	result[a_Length] = '\0';
}

// Generate a float between a_Min and a_Max
void generateFloat( std::size_t a_Length, char* result )
{
	int decimalPointIdx = randInt( 0, static_cast<int>( a_Length ) - 2 );
	index_t startIdx = 0;

	if ( decimalPointIdx == 0 )
	{
		result[0] = '0';
		result[1] = '.';
		startIdx = 2;
	}

	for ( index_t i = startIdx; i < a_Length + startIdx; i++ )
	{
		if ( decimalPointIdx == i )
			result[i] = '.';
		else
			result[i] = static_cast<char>( randInt( 0, 9 ) ) + '0';	// Get random number between 0 and 10
	}
	result[a_Length] = '\0';
}

Athena::Number generateRandomNumber( std::size_t a_Precision, bool isFloating )
{
	std::size_t numDigits = (a_Precision * 3) / 10 + 1;	// Estimate number of digits needed for given precision (log10(2) ≈ 0.30103)
	std::vector<char> buffer( numDigits + 2 );	// +1 for possible decimal point, +1 for null terminator
	if ( isFloating )
		generateFloat( numDigits, buffer.data() );
	else
		generateInt( numDigits, buffer.data() );
	return Athena::Number( buffer.data(), a_Precision );
}