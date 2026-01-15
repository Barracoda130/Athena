#pragma once
#include <vector>
#include <string>

/*
 * Numbers stored using number_str_t are stored backwards
 * ie. the number 1345 is stored as [5,4,3,1]
 * This is to make the mathmatic calculations easier as the size of the result will
 * often need to change at the largest digit, so that digit is now on the RHS.
 * Negative signs are still placed at the front of the number however.
 * 
 * Strings have been implemented using std:vector rather than std:string
 * because of significant performance improvement. The conversion time
 * to convert from string to vector is worth it considering the number
 * of calculations required to be performed.
 */

enum sign_t { POSITIVE, NEGATIVE, NOT_A_NUMBER };		// The sign can be used to determine whether is a valid number or not

struct number_str_t
{
	std::vector<char> data;
	sign_t sign;
};
// typedef std::vector<char> number_str_t;

number_str_t toNumberStr( const std::string& a_Str );
std::string fromNumberStr( const number_str_t& a_Num );

number_str_t strAdd( const number_str_t& a_N1, const number_str_t& a_N2 );