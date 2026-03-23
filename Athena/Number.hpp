#pragma once
#include <vector>
#include <string>
#include "MpfrInclude.hpp"

namespace Athena
{
	typedef mpfr_prec_t precision_t;
	typedef mpfr_rnd_t round_t;
	typedef mpfr_exp_t exponent_t;
	typedef mp_limb_t limb_t;
	typedef mpfr_sign_t sign_t;
	

	class Number
	{
	public:
		Number();
		Number( precision_t a_Precision );
		Number( const std::string& a_Value, precision_t a_Precision );
		Number( const long long a_Value, precision_t a_Precison );
		Number( const Number& a_Value, round_t a_Round ) : Number( a_Value.m_Value, a_Round ) {};
		Number( const mpfr_t a_Value, round_t a_Round );

		// Destructor
		~Number();

		// Set methods will not adjust the precision of this number
		void set( const Number& a_Value, round_t a_Round );
		void set( const mpfr_t a_Value, round_t a_Round );
		void set( const std::string& a_Value, round_t a_Round );
		void set( const long long a_Value, round_t a_Round );

		std::string str() const;

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
       Number& operator=( const Number& a_Other );

		// Friend methods
		friend void add( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );
		friend void sub( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );
		friend void mul( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );
		friend void div( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );

		// Other
		friend std::ostream& operator<<( std::ostream& os, const Number& a_Num );

	private:
		exponent_t getExp() const { return m_Value->_mpfr_exp; };
		precision_t getPrec() const { return m_Value->_mpfr_prec; };
		sign_t getSign() const { return m_Value->_mpfr_sign; };
		const limb_t *getLimbs() const { return m_Value->_mpfr_d; };
		limb_t* getLimbs() { return m_Value->_mpfr_d; };

		mpfr_t m_Value;
	};

	// Namespace-scope declarations so callers can use Athena::add / sub / mul / div
	void add( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );
	void sub( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );
	void mul( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );
	void div( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round );

	// Benchmark function for comparing add_n implementations
	void benchmark_add_n_implementations( size_t limbCount = 100, size_t iterations = 10000 );

	// Namespace-scope declaration for the stream operator
	std::ostream& operator<<( std::ostream& os, const Number& a_Num );
}