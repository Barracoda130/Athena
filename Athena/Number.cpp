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
		mpfr_init( m_Value );
	}

	Number::Number( precision_t a_Precision )
	{
		// Maybe use macro in future: MPFR_DECL_INIT
		mpfr_init2( m_Value, a_Precision );
	}

	Number::Number( const std::string& a_Value, precision_t a_Precision )
	{
		mpfr_init2( m_Value, a_Precision );
		set( a_Value, MPFR_RNDN );
	}

	Number::Number( const long long a_Value, precision_t a_Precision )
	{
		setPrecision( a_Precision );
		set( a_Value );
	}

	// Set methods
	void Number::set( const Number& a_Value )
	{
		mpfr_set( m_Value, a_Value.m_Value, MPFR_RNDN );
	}

	void Number::set( const std::string& a_Value, round_t a_Round )
	{
		mpfr_set_str( m_Value, a_Value.c_str(), 10, a_Round );
	}

	void Number::set( const long long a_Value )
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
		set( a_Other );
		return *this;
	}

// Private

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
