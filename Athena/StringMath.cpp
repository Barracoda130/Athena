#include "StringMath.hpp"

#include "Common.hpp"

number_str_t toNumberStr( const std::string& a_Str )
{
	number_str_t result;

	index_t minIndex = 0;
	// If the first index of both numbers is "-" then treat min index as 1
	if ( a_Str[0] == '-' )
	{
		result.sign = NEGATIVE;
		minIndex = 1;
	}
	else
		result.sign = POSITIVE;

	result.data.resize( a_Str.size() - minIndex );	// - min index to account for possible '-'

	for ( index_t i = minIndex; i < a_Str.size(); i++ )
	{
		result.data[result.data.size() + minIndex - i - 1] = a_Str[i];
	}
	return result;
}

std::string fromNumberStr( const number_str_t& a_Num )
{
	std::string result;
	
	index_t minIndex = 0;
	if ( a_Num.sign == NEGATIVE )
	{
		result.resize( a_Num.data.size() + 1 );	// +1 to account for '-'
		result[0] = '-';
		minIndex = 1;
	}
	else
	{
		result.resize( a_Num.data.size() );
	}

	for ( index_t i = minIndex; i < a_Num.data.size() + minIndex; i++ )
	{
		result[result.size() + minIndex - i - 1] = a_Num.data[i - minIndex];
	}
	return result;
}

number_str_t strAdd( const number_str_t& a_N1, const number_str_t& a_N2 )
{
	// If one number is negative and the other positive, then we need to do subtraction instead


	const char zero = '0';

	number_str_t result;
	result.sign = a_N1.sign;	// Both N1 and N2 have the same sign at this point
	result.data.resize( std::max( a_N1.data.size(), a_N2.data.size() ) );

	char carry = 0;
	index_t i;
	for ( i = 0; i < a_N1.data.size() && i < a_N2.data.size(); i++ )
	{
		char current = a_N1.data[i] + a_N2.data[i] + carry - zero;
		if ( current > zero + 9 )
		{
			current -= 10;
			carry = 1;
		}
		else
		{
			carry = 0;
		}
		result.data[i] = current;
	}

	if ( i < a_N1.data.size() )
	{
		while ( i < a_N1.data.size() && carry != 0 )
		{
			char current = a_N1.data[i] + carry;
			if ( current > zero + 9 )
			{
				current -= 10;
				carry = 1;
			}
			else
			{w
				carry = 0;
			}
			result.data[i] = current;
		}

		while ( i < a_N1.data.size() )
			result.data[i] = a_N1.data[i];
	}

	else if ( i < a_N2.data.size() )
	{
		while ( i < a_N2.data.size() && carry != 0 )
		{
			char current = a_N2.data[i] + carry;
			if ( current > zero + 9 )
			{
				current -= 10;
				carry = 1;
			}
			else
			{
				carry = 0;
			}
			result.data[i] = current;
		}

		while ( i < a_N2.data.size() )
			result.data[i] = a_N2.data[i];
	}

	// If we got to this stage and the carry is not 0, then we need to increase the length of the result by 1
	if ( carry != 0 )
		result.data.emplace_back( zero + 1 );

	return result;
}

number_str_t strSub( const number_str_t& a_N1, const number_str_t& a_N2 )
{
	// This will subtract N2 from N1
	number_str_t result;
	result.data.resize( std::max( a_N1.data.size(), a_N2.data.size() ) );


	return result;
}