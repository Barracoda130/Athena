#pragma once
#include <cstddef>
#include <fstream>
#include <string>

#include "common.hpp"

namespace TestData
{
	// Types
	typedef enum { POSITIVE, NEGATIVE, RANDOM } sign_t;
	struct column_data_t
	{
		sign_t sign;
		bool isFloating;
		char name[8];
	};

	const char* getFileName();
	void generateData(std::size_t a_NumRows);

	class FileReader
	{
	public:
		FileReader( bool a_AFloating, sign_t a_ASign,
					bool a_BFloating, sign_t a_BSign );

		~FileReader();

		void reset( bool a_AFloating, sign_t a_ASign,
					bool a_BFloating, sign_t a_BSign );

		const char* getNextA();
		const char* getNextB();

		bool eofReached() { return m_EOFReached; };

	private:
		void readNextLine();

		std::ifstream m_File;
		index_t m_ColAIdx;
		index_t m_ColBIdx;
		index_t m_CurrentRowIdx;

		bool m_ABeenRead;
		bool m_BBeenRead;
		std::string m_AData;
		std::string m_BData;

		bool m_EOFReached;

	};
}

