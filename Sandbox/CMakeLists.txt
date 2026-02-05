# CMakeList.txt : CMake project for Sandbox, include source and define
# project specific logic here.
#

pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

# Add source to this project's executable.
add_executable (Sandbox "Sandbox.cpp" "Sandbox.h")

target_link_libraries(Sandbox PRIVATE PkgConfig::mpfr)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Sandbox PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
