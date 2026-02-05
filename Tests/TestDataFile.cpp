#include "TestDataFile.hpp"

#include <array>
#include <string>
#include <cstdlib>
#include <cassert>
#include <filesystem>

#define MIN_NUM_LENGTH 2
#define MAX_NUM_LENGTH 10

namespace TestData
{
	// Static vars
	const char s_FileName[] = "TestData.csv";
	std::array<column_data_t, 12> s_ColumnHeaders = { {
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
		return s_FileName;
	}

	void generateInt( std::size_t a_Length, char* result )
	{
		for ( index_t i = 0; i < a_Length; i++ )
		{
			result[i] = static_cast< char >( randInt( 0, 9 ) ) + '0';	// Get random number between 0 and 10
		}
		result[a_Length] = '\0';
	}

	// Generate a float between a_Min and a_Max
	void generateFloat( std::size_t a_Length, char* result )
	{
		int decimalPointIdx = randInt( 0, static_cast<int>(a_Length) - 2 );
		index_t startIdx = 0;

		if ( decimalPointIdx == 0 )
		{
			result[0] = '0';
			result[1] = '.';
			startIdx = 2;
		}

		for ( index_t i = startIdx; i < a_Length + startIdx; i++ )
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
		std::ofstream f( s_FileName );

		// Add all headers to the file
		for ( index_t i = 0; i < s_ColumnHeaders.size(); i++ )
		{
			f << s_ColumnHeaders[i].name;

			if ( i == s_ColumnHeaders.size() - 1 )
				f << std::endl;
			else
				f << ",";
		}
		

		// Add numbers for columns
		for ( int i = 0; i < a_NumRows; i++ )
		{
			for ( index_t columnIdx = 0; columnIdx < s_ColumnHeaders.size(); columnIdx++ )
			{
				char buf[MAX_NUM_LENGTH + 3];

				std::size_t numLength = static_cast< std::size_t >( randInt( MIN_NUM_LENGTH, MAX_NUM_LENGTH ) );
				if ( !s_ColumnHeaders[columnIdx].isFloating )
					generateInt( numLength, buf );
				else
					generateFloat( numLength, buf );;// generate a float

				// If negative, write the negative sign to the file
				if ( s_ColumnHeaders[columnIdx].sign == NEGATIVE )
					f << "-";
				else if ( s_ColumnHeaders[columnIdx].sign == RANDOM )
					if ( std::rand() % 2 == 0 )		// If sign is random, then 50% chance of it being negative
						f << "-";

				f << buf;
				if ( columnIdx == s_ColumnHeaders.size() - 1 )
					f << std::endl;
				else
					f << ",";
			}
		}
		f.close();
	}

	std::string getNthColumn( std::string& a_Input, index_t n )
	{
		// If looking for first column just return that
		if ( n == 0 )
			return a_Input.substr( 0, a_Input.find( "," ) );

		std::string result;
		index_t currColIdx = 0;

		index_t start = 0;
		index_t end = 0;

		for ( index_t i = 0; i < a_Input.size(); i++ )
		{
			if ( a_Input[i] == ',' )
			{
				currColIdx++;
				
				if ( currColIdx == n )
					start = i;
				else if ( currColIdx == n + 1 )
					end = i;
			}
		}

		return a_Input.substr( start + 1, end - start ); // +1 to remove the ,
	}


	
	FileReader::FileReader( bool a_AFloating, sign_t a_ASign, bool a_BFloating, sign_t a_BSign )
	{
		reset( a_AFloating, a_ASign, a_BFloating, a_BSign );
	}

	FileReader::~FileReader()
	{
		m_File.close();
	}

	void FileReader::reset( bool a_AFloating, sign_t a_ASign, bool a_BFloating, sign_t a_BSign )
	{
		m_ABeenRead = false;
		m_BBeenRead = false;
		m_ColAIdx = 0;
		m_ColBIdx = 0;

		for ( index_t i = 0; i < s_ColumnHeaders.size(); i++ )
		{
			column_data_t currCol = s_ColumnHeaders[i];
			if ( !m_ABeenRead &&				// Not assigned yet
				currCol.isFloating == a_AFloating &&
				currCol.sign == a_ASign )
			{
				m_ColAIdx = i;
				m_ABeenRead = true;
			}
			else if ( !m_BBeenRead &&			// Not assigned yet
				currCol.isFloating == a_BFloating &&
				currCol.sign == a_BSign )
			{
				m_ColBIdx = i;
				m_BBeenRead = true;
			}
		}

		assert( m_ColAIdx != m_ColBIdx &&	// They shouldn't be pointing to the same column
			m_ColAIdx >= 0 &&			// They each should have a column
			m_ColBIdx >= 0 );

		m_EOFReached = false;

		m_File = std::ifstream( s_FileName );
		m_CurrentRowIdx = 0;

		readNextLine();		// Skip the header
		readNextLine();
	}

	void FileReader::readNextLine()
	{
		if ( m_EOFReached )
			return;

		std::string buf;
		if ( !getline( m_File, buf ) )
		{
			m_EOFReached = true;
			m_File.close();
			m_AData = "";
			m_BData = "";
			return;
		}

		m_AData = getNthColumn( buf, m_ColAIdx );
		m_BData = getNthColumn( buf, m_ColBIdx );

		m_ABeenRead = false;
		m_BBeenRead = false;
	}

	const char* FileReader::getNextA()
	{
		// If This value has already been read get a new one
		if ( m_ABeenRead && !m_EOFReached )
			readNextLine();

		m_ABeenRead = true;
		return m_AData.c_str();
	}

	const char* FileReader::getNextB()
	{
		// If This value has already been read get a new one
		if ( m_BBeenRead && !m_EOFReached )
			readNextLine();

		m_BBeenRead = true;
		return m_BData.c_str();
	}
}



