pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

include_directories("../Common")

add_library(Athena)

target_sources(Athena
	PRIVATE
		"Number.hpp"
		"Number.cpp"
		"Athena.cpp"
		"NumberAdd.cpp"
	PUBLIC
		FILE_SET HEADERS
		FILES
			"Athena.hpp"

)

target_include_directories(Athena PRIVATE ${CMAKE_SOURCE_DIR}/Common)

target_link_libraries(Athena PRIVATE PkgConfig::mpfr)