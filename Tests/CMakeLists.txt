pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

include_directories("../Athena")
include_directories("../Common")

# Add source to this project's executable.
add_executable (Tests "RunTests.cpp" "MpfrInclude.hpp" "Types.hpp" "TestDataFile.cpp" "TestDataFile.hpp")

target_link_libraries(Tests PRIVATE PkgConfig::mpfr)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Tests PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
