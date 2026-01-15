#pragma once
#include <vector>
#include <string>
#include "MpfrInclude.hpp"
#include "StringMath.hpp"

namespace Athena
{
	typedef long exponent_t;
	typedef std::size_t precision_t;
	typedef unsigned long long mantissa_t;

	// Roundings are taken from MPFR
	enum round_t
	{
		RNDD,	// Round down towards -infinity
		RNDN,	// Round to the nearest value (up or down)
		RNDU,	// Round up towards +infinity
		RNDZ,	// Round towards zero (truncate)
		RNDA	// Round away from zero
	};

	class Number
	{
	public:
		Number();
		Number( precision_t a_Precision );
		Number( const std::string& a_Value, precision_t a_Precision );
		Number( const long long a_Value, precision_t a_Precison );

		// Set methods will not adjust the precision of this number
		void set( const Number& a_Value );
		void set( const std::string& a_Value, round_t a_Round );
		void set( const long long a_Value );

		bool beenInitialised() const { return m_Precision > 0; };

		//const std::vector<mantissa_t>& getMantissa() const { return m_Mantissa; };

		// Operator overloads
		// Comparison
		bool operator== ( const mpfr_t& a_Other ) const;	// Used for testing
		bool operator== ( const Number& a_Other ) const;
		bool operator!= ( const Number& a_Other ) const;
		bool operator> ( const Number& a_Other ) const;
		bool operator< ( const Number& a_Other ) const;
		bool operator>= ( const Number& a_Other ) const;
		bool operator<= ( const Number& a_Other ) const;

		// Assignment
		Number& operator=( Number& a_Other );

		// Friend methods
		friend void add( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );
		friend void sub( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );
		friend void mult( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );
		friend void div( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );

	private:
		void setPrecision( precision_t a_Precision );

		exponent_t m_Exp;
		sign_t m_Sign;
		precision_t m_Precision;
		std::vector<mantissa_t> m_Mantissa;
	};

	
}