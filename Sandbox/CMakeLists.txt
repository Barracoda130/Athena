pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

include_directories("../Athena")
include_directories("../Common")

# Add source to this project's executable.
add_executable (Sandbox "main.cpp" "MpfrInclude.hpp" "MpfrPlayground.cpp" "MpfrPlayground.hpp")

target_link_libraries(Sandbox PRIVATE PkgConfig::mpfr)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Sandbox PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
