#include "Number.hpp"

#include <cassert>
#include <format>
#include <iostream>
#include "Common.hpp"

namespace Athena
{
	Number::Number()
	{
		mpfr_init( m_Value );
	}

	Number::Number( precision_t a_Precision )
	{
		// Maybe use macro in future: MPFR_DECL_INIT
		// Will need to modify the destructor if i do this
		mpfr_init2( m_Value, a_Precision );
	}

	Number::Number( const std::string& a_Value, precision_t a_Precision )
	{
		mpfr_init2( m_Value, a_Precision );
		set( a_Value, MPFR_RNDN );
	}

	Number::Number( const long long a_Value, precision_t a_Precision )
	{
		mpfr_init2( m_Value, a_Precision );
		set( a_Value, MPFR_RNDN );
	}

	Number::Number( const mpfr_t a_Value, round_t a_Round )
	{
		mpfr_init2( m_Value, mpfr_get_prec( a_Value ) );
		set( a_Value, a_Round );
	}

	Number::~Number()
	{
		mpfr_clear( m_Value );
	}

	// Set methods
	void Number::set( const Number& a_Value, round_t a_Round )
	{
		set( a_Value.m_Value, a_Round );
	}

	void Number::set( const mpfr_t a_Value, round_t a_Round )
	{
		mpfr_set( m_Value, a_Value, a_Round );
	}

	void Number::set( const std::string& a_Value, round_t a_Round )
	{
		mpfr_set_str( m_Value, a_Value.c_str(), 10, a_Round );
	}

	void Number::set( const long long a_Value, round_t a_Round )
	{
		mpfr_set_uj( m_Value, static_cast<unsigned long>( a_Value ), a_Round );
	}

	bool Number::operator==( const Number& a_Other ) const
	{
		return static_cast<bool>( mpfr_equal_p( m_Value, a_Other.m_Value ) );
	}

	// For testing
	bool Number::operator==( const mpfr_t& a_Other ) const
	{
		bool result = static_cast<bool>( mpfr_equal_p( m_Value, a_Other ) );
		if ( !result )
		{
			Number tmp( a_Other, MPFR_RNDN );
			std::cout << *this << " != " << tmp << std::endl;
		}
		return result;
	}

	bool Number::operator!=( const Number& a_Other ) const
	{
		return !static_cast<bool>( mpfr_equal_p( m_Value, a_Other.m_Value ) );
	}

	bool Number::operator>( const Number& a_Other ) const
	{
		return static_cast<bool>( mpfr_greater_p( m_Value, a_Other.m_Value ) );
	}

	bool Number::operator<( const Number& a_Other ) const
	{
		return static_cast<bool>( mpfr_greater_p( a_Other.m_Value, m_Value ) );
	}

	bool Number::operator>=( const Number& a_Other ) const
	{
		return static_cast<bool>( mpfr_greaterequal_p( m_Value, a_Other.m_Value ) );
	}

	bool Number::operator<=( const Number& a_Other ) const
	{
		return static_cast<bool>( mpfr_greaterequal_p( a_Other.m_Value, m_Value ) );
	}

	Number& Number::operator=( Number& a_Other )
	{
		set( a_Other, MPFR_RNDN );
		return *this;
	}

// Private

// Friend 
	std::ostream& operator<<( std::ostream& os, const Number& a_Num )
	{
		mpfr_exp_t exponent;
		char* mantissa = mpfr_get_str( nullptr, &exponent, 10, 0, a_Num.m_Value, MPFR_RNDN );
		char initial = mantissa[0];
		mantissa++;

		return os << std::format( "{}.{}e{}", initial, mantissa, exponent - 1);
	}

	// Non member methods

	void sub( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round )
	{
		mpfr_sub( a_Result.m_Value, a_Num1.m_Value, a_Num2.m_Value, a_Round );
	}

	void mul( Number & a_Result, const Number & a_Num1, const Number & a_Num2, round_t a_Round )
	{
		mpfr_mul( a_Result.m_Value, a_Num1.m_Value, a_Num2.m_Value, a_Round );
	}

	void div( Number & a_Result, const Number & a_Num1, const Number & a_Num2, round_t a_Round )
	{
		mpfr_div( a_Result.m_Value, a_Num1.m_Value, a_Num2.m_Value, a_Round );
	}

}
