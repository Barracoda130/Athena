pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
find_package(benchmark CONFIG QUIET)

# Add source to this project's executable.
add_executable (Benchmarking)



target_sources(Benchmarking
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/GoogleBenchmark.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/GoogleBenchmark.hpp"
        ${COMMON_SRC}
        ${COMMON_HDR})

# Add the Common directory to the target include paths so headers can be referenced without "../"
target_include_directories(Benchmarking PRIVATE ${CMAKE_SOURCE_DIR}/Common)

target_link_libraries(Benchmarking PRIVATE PkgConfig::mpfr)
target_link_libraries(Benchmarking PRIVATE Athena)


target_link_libraries(Benchmarking PRIVATE benchmark::benchmark)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Benchmarking PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
