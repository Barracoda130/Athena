#include "TestDataFile.hpp"

#include <array>
#include <fstream>
#include <string>
#include <cstdlib>

#define MIN_NUM_LENGTH 10
#define MAX_NUM_LENGTH 40

namespace TestData
{
	typedef enum {POSITIVE, NEGATIVE, RANDOM} sign_t;
	struct column_data_t
	{
		sign_t sign;
		bool isFloating;
		char name[8];
	};

	const char fileName[] = "TestData.csv";
	std::array<column_data_t, 12> columnHeaders = { {
		{POSITIVE, false, "posInt1"},
		{POSITIVE, false, "posInt2"},
		{NEGATIVE, false, "negInt1"},
		{NEGATIVE, false, "negInt2"},
		{RANDOM, false, "rndInt1"},
		{RANDOM, false, "rndInt2"},
		{POSITIVE, true, "posFlt1"},
		{POSITIVE, true, "posFlt2"},
		{NEGATIVE, true, "negFlt1"},
		{NEGATIVE, true, "negFlt2"},
		{RANDOM, true, "rndFlt1"},
		{RANDOM, true, "rndFlt2"},
	} };

	// Generate a random int between a_Low and a_High INCLUSIVE
	int randInt( int a_Low, int a_High )
	{
		return ( std::rand() % ( a_High + 1 - a_Low ) ) + a_Low;
	}

	const char* getFileName()
	{
		return fileName;
	}

	void generateInt( std::size_t a_Length, char* result )
	{
		for ( std::size_t i = 0; i < a_Length; i++ )
		{
			result[i] = static_cast< char >( randInt( 0, 9 ) ) + '0';	// Get random number between 0 and 10
		}
		result[a_Length] = '\0';
	}

	// Generate a float between a_Min and a_Max
	void generateFloat( std::size_t a_Length, char* result )
	{
		int decimalPointIdx = randInt( 0, static_cast<int>(a_Length) - 2 );
		std::size_t startIdx = 0;

		if ( decimalPointIdx == 0 )
		{
			result[0] = '0';
			result[1] = '.';
			startIdx = 2;
		}

		for ( std::size_t i = startIdx; i < a_Length + startIdx; i++ )
		{
			if ( decimalPointIdx == i )
				result[i] = '.';
			else
				result[i] = static_cast< char >( randInt( 0, 9 ) ) + '0';	// Get random number between 0 and 10
		}
		result[a_Length] = '\0';
	}

	// Generate aNumRows rows of test data
	void generateData( std::size_t a_NumRows )
	{
		std::ofstream f( fileName );

		// Add all headers to the file
		for ( std::size_t i = 0; i < columnHeaders.size(); i++ )
		{
			f << columnHeaders[i].name;

			if ( i == columnHeaders.size() - 1 )
				f << std::endl;
			else
				f << ",";
		}
		

		// Add numbers for columns
		for ( int i = 0; i < a_NumRows; i++ )
		{
			for ( std::size_t columnIdx = 0; columnIdx < columnHeaders.size(); columnIdx++ )
			{
				char buf[MAX_NUM_LENGTH + 3];

				std::size_t numLength = static_cast< std::size_t >( randInt( MIN_NUM_LENGTH, MAX_NUM_LENGTH ) );
				if ( !columnHeaders[columnIdx].isFloating )
					generateInt( numLength, buf );
				else
					generateFloat( numLength, buf );;// generate a float

				// If negative, write the negative sign to the file
				if ( columnHeaders[columnIdx].sign == NEGATIVE )
					f << "-";
				else if ( columnHeaders[columnIdx].sign == RANDOM )
					if ( std::rand() % 2 == 0 )		// If sign is random, then 50% chance of it being negative
						f << "-";

				f << buf;
				if ( columnIdx == columnHeaders.size() - 1 )
					f << std::endl;
				else
					f << ",";
			}
		}


		f.close();
	}

	
}



