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

	Number::Number( const Number& a_Value, round_t a_Round )
	{
		mpfr_init2( m_Value, mpfr_get_prec( a_Value.m_Value ) );
		set( a_Value, a_Round );
	}

	Number::~Number()
	{
		mpfr_clear( m_Value );
	}

	// Set methods
	void Number::set( const Number& a_Value, round_t a_Round )
	{
		mpfr_set( m_Value, a_Value.m_Value, MPFR_RNDN );
	}

	void Number::set( const std::string& a_Value, round_t a_Round )
	{
		mpfr_set_str( m_Value, a_Value.c_str(), 10, a_Round );
	}

	void Number::set( const long long a_Value, round_t a_Round )
	{
		mpfr_set_uj( m_Value, static_cast<unsigned long>( a_Value ), MPFR_RNDN );
	}

	bool Number::operator==( const Number& a_Other ) const
	{
		return false;
	}

	// For testing
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
		set( a_Other, MPFR_RNDN );
		return *this;
	}

// Private

	// Non member methods
	void add( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round )
	{
		// This line should be removed for efficiency later
		assert( true );

		//
	}

// Friend 
	std::ostream& operator<<( std::ostream& os, Number& a_Num )
	{
		mpfr_exp_t exponent;
		char* mantissa = mpfr_get_str( nullptr, &exponent, 10, 0, a_Num.m_Value, MPFR_RNDN );
		char initial = mantissa[0];
		mantissa++;

		return os << std::format( "{}.{}e{}", initial, mantissa, exponent - 1);
	}

	void sub( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round )
	{}

	void mult( Number & a_Result, const Number & a_Num1, const Number & a_Num2, round_t a_Round )
	{}

	void div( Number & a_Result, const Number & a_Num1, const Number & a_Num2, round_t a_Round )
	{}

}
