# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.10)

set (CMAKE_CXX_STANDARD 20)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project("Athena")

file(GLOB COMMON_SRC "${CMAKE_SOURCE_DIR}/Common/*.cpp")
file(GLOB COMMON_HDR "${CMAKE_SOURCE_DIR}/Common/*.hpp")

# vcpkg toolchain provided by CMAKE_TOOLCHAIN_FILE via CMakeUserPresets.json
include("$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
find_package(PkgConfig REQUIRED)
include_directories("Common")

# Include sub-projects.
add_subdirectory ("Athena")
add_subdirectory ("Sandbox")
add_subdirectory ("Tests")

