pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

# Add source to this project's executable.
add_executable (Sandbox)



target_sources(Sandbox
    PRIVATE
        "main.cpp"
        "MpfrPlayground.cpp"
        "MpfrPlayground.hpp"
        ${COMMON_SRC}
        ${COMMON_HDR}
 "CodeProfiling.cpp")

# Add the Common directory to the target include paths so headers can be referenced without "../"
target_include_directories(Sandbox PRIVATE ${CMAKE_SOURCE_DIR}/Common)

target_link_libraries(Sandbox PRIVATE PkgConfig::mpfr)
target_link_libraries(Sandbox PRIVATE Athena)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Sandbox PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
