pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

include_directories("../Common")
include_directories("mpfr-required")

add_library(Athena)

# Collect mpfr implementation sources and headers from the bundled directory so
# individual relative paths aren't needed in target_sources.
file(GLOB MPFR_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/mpfr-required/*.c")
file(GLOB MPFR_HDRS "${CMAKE_CURRENT_SOURCE_DIR}/mpfr-required/*.h")

target_sources(Athena
    PRIVATE
        "Number.hpp"
        "Number.cpp"
        "Athena.cpp"
        "NumberAdd.cpp"
        ${MPFR_SRCS}
        ${MPFR_HDRS}
        ${COMMON_SRC}
        ${COMMON_HDR}
    PUBLIC
        FILE_SET HEADERS
        FILES
            "Athena.hpp"

   "RShift.cpp")

target_include_directories(Athena PRIVATE ${CMAKE_SOURCE_DIR}/Common)

# Add mpfr-required to include directories so headers inside it are found by
# the compiler when using #include "mpfr-*.h" without specifying a path.
target_include_directories(Athena PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/mpfr-required")

target_link_libraries(Athena PRIVATE PkgConfig::mpfr)
