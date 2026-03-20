pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

# Add source to this project's executable.
add_executable (Benchmarking)



target_sources(Benchmarking
    PRIVATE
        "main.cpp"
        "Benchmark.cpp"
        "Benchmark.hpp"
        "BenchmarkFunction.cpp"
        "BenchmarkFunction.hpp"
        ${COMMON_SRC}
        ${COMMON_HDR} "BenchmarkFunction.cpp")

# Add the Common directory to the target include paths so headers can be referenced without "../"
target_include_directories(Benchmarking PRIVATE ${CMAKE_SOURCE_DIR}/Common)

target_link_libraries(Benchmarking PRIVATE PkgConfig::mpfr)
target_link_libraries(Benchmarking PRIVATE Athena)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Benchmarking PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
