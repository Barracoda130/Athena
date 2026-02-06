pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

# Add source to this project's executable.
add_executable(Tests)

target_sources(Tests
	PRIVATE
		"RunTests.cpp"
		"TestDataFile.cpp"
		"TestDataFile.hpp"
)

target_include_directories(Tests PRIVATE ${CMAKE_SOURCE_DIR}/Common)

target_link_libraries(Tests PRIVATE PkgConfig::mpfr)
target_link_libraries(Tests PRIVATE Athena)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Tests PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
