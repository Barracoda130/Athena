pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

# Use Catch2's CMake config to get the target that supplies a main implementation
find_package(Catch2 CONFIG REQUIRED)

# Add source to this project's executable.
add_executable(Tests)
add_executable(Catch2Sandbox)

target_sources(Tests
	PRIVATE
		"RunTests.cpp"
		"TestDataFile.cpp"
		"TestDataFile.hpp"
		${COMMON_SRC}
        ${COMMON_HDR} )

target_sources(Catch2Sandbox
	PRIVATE
		"RunCatch2.cpp")

target_include_directories(Tests PRIVATE ${CMAKE_SOURCE_DIR}/Common)

target_link_libraries(Tests PRIVATE PkgConfig::mpfr)
target_link_libraries(Tests PRIVATE Catch2::Catch2WithMain)
target_link_libraries(Tests PRIVATE Athena)

target_link_libraries(Catch2Sandbox PRIVATE Catch2::Catch2WithMain)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Tests PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
