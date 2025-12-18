#include "TestDataFile.hpp"

#include <array>
#include <fstream>

static const char s_FileName[] = "TestData.csv";
static std::array<const char[7], 12> s_ColumnHeaders = {
	"posInt", "posInt",
	"negInt", "negInt",
	"rndInt", "rndInt",
	"posFlt", "posFlt",
	"negFlt", "negFlt",
	"rndFlt", "rndFlt"
};

namespace TestData
{
	const char* getFileName()
	{
		return s_FileName;
	}

	void generateData( std::size_t a_NumRows )
	{
		// Generate aNumRows rows of test data

		// Add all headers to the file
		std::ofstream f( s_FileName );
		
		for ( int i = 0; i < s_ColumnHeaders.size(); i++ )
		{
			f << s_ColumnHeaders[i];
		}
		f << std::endl;
		f.close();
	}
}



