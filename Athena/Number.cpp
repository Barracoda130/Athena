#include "Number.hpp"

#include <cassert>
#include "Common.hpp"

namespace
{
	constexpr std::size_t mantissa_tBase10Length = 0;
}

namespace Athena
{
	Number::Number()
	{
		m_Precision = 0;
		m_Exp = 0;
		m_Sign = NOT_A_NUMBER;
		
	}

	Number::Number( precision_t a_Precision )
	{
		m_Exp = 0;
		setPrecision( a_Precision );
		m_Sign = NOT_A_NUMBER;
	}

	Number::Number( const std::string& a_Value, precision_t a_Precision )
	{
		setPrecision( a_Precision );
		set( a_Value, RNDD );
	}

	Number::Number( const long long a_Value, precision_t a_Precision )
	{
		setPrecision( a_Precision );
		set( a_Value );
	}

	// Set methods
	void Number::set( const Number& a_Value )
	{
		m_Exp = a_Value.m_Exp;
		m_Sign = a_Value.m_Sign;
		m_Mantissa = a_Value.m_Mantissa;
	}

	void Number::set( const std::string& a_Value, round_t a_Round )
	{
		// Temporary version!!
		long long llValue = std::stoll( a_Value );
		m_Exp = 0;
		m_Sign = llValue >= 0 ? POSITIVE : NEGATIVE;
		m_Mantissa[0] = static_cast<mantissa_t>(llValue);

		// Initial slow version
		m_Exp = static_cast<exponent_t>( a_Value.size() );
		m_Sign = a_Value[0] == '-' ? NEGATIVE : POSITIVE;

		// If the number is too large for the prevision we are allowed then remove the least significant bits
		// We will fill right to left, but if the input number is larger than we can represent
		// Then we will trim the string
		for ( int pow = 0; ; pow++ )
		{
			index_t i = 0;
		}
	}

	void Number::set( const long long a_Value )
	{
		m_Exp = 0;
		m_Sign = a_Value >= 0 ? POSITIVE : NEGATIVE;
		m_Mantissa[0] = static_cast<mantissa_t>( a_Value );
	}

	bool Number::operator==( const Number& a_Other ) const
	{
		return false;
	}

	bool Number::operator==( const mpfr_t& a_Other ) const
	{
		return false;
	}

	bool Number::operator!=( const Number& a_Other ) const
	{
		return false;
	}

	bool Number::operator>( const Number& a_Other ) const
	{
		return false;
	}

	bool Number::operator<( const Number& a_Other ) const
	{
		return false;
	}

	bool Number::operator>=( const Number& a_Other ) const
	{
		return false;
	}

	bool Number::operator<=( const Number& a_Other ) const
	{
		return false;
	}

	Number& Number::operator=( Number& a_Other )
	{
		set( a_Other );
		return *this;
	}

// Private
	void Number::setPrecision( precision_t a_Precision )
	{
		// Precision is given in bits
		// Will always set precision to the nearest 
		m_Precision = a_Precision;
		std::size_t mantissaSizeBits = sizeof( mantissa_t ) * 8;


		m_Mantissa.resize( a_Precision / mantissaSizeBits + 1 );
	}

	// Non member methods
	void add( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round )
	{
		// This line should be removed for efficiency later
		assert( a_Result.beenInitialised() && a_Num1.beenInitialised() && a_Num2.beenInitialised() );

		//
	}

	void sub( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round )
	{}

	void mult( Number & a_Result, const Number & a_Num1, const Number & a_Num2, round_t a_Round )
	{}

	void div( Number & a_Result, const Number & a_Num1, const Number & a_Num2, round_t a_Round )
	{}

}
