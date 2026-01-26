#pragma once
#include <vector>
#include <string>
#include "MpfrInclude.hpp"

namespace Athena
{
	typedef mpfr_prec_t precision_t;
	typedef mpfr_rnd_t round_t;

	class Number
	{
	public:
		Number();
		Number( precision_t a_Precision );
		Number( const std::string& a_Value, precision_t a_Precision );
		Number( const long long a_Value, precision_t a_Precison );
		Number( const Number& a_Value, round_t a_Round );

		// Destructor
		~Number();

		// Set methods will not adjust the precision of this number
		void set( const Number& a_Value, round_t a_Round );
		void set( const std::string& a_Value, round_t a_Round );
		void set( const long long a_Value, round_t a_Round );

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

		// Other
		friend std::ostream& operator<<( std::ostream& os, Number& a_Num );

	private:

		mpfr_t m_Value;
	};

	
}