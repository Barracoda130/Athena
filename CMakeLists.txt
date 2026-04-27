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

option(ATHENA_ENABLE_FAST_MATH "Enable aggressive floating-point optimizations" OFF)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Compiler-specific optimization flags
if(MSVC)
    # MSVC optimization flags
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /Ob2 /Oi /Ot /GL")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2 /Ob2 /Oi /Ot /GL")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
    set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS_RELEASE} /LTCG")

    # Additional MSVC optimizations for modern x64 CPUs (e.g. Zen 4)
    add_compile_options(
        $<$<CONFIG:Release>:/arch:AVX512>
        $<$<CONFIG:Release>:/favor:AMD64>
    )

    if(ATHENA_ENABLE_FAST_MATH)
        add_compile_options($<$<CONFIG:Release>:/fp:fast>)
    endif()
else()
    # GCC/Clang optimization flags
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ flags (Release): ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "C flags (Release): ${CMAKE_C_FLAGS_RELEASE}")

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
add_subdirectory ("Benchmarking")

