# examples/example/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(CppEOSDemo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# This example is part of the main project build.
# Link directly to the targets defined in the parent project.
# Ensure these targets (EOSCore, EquationOfStateV1_lib) are defined before this point.
if(NOT TARGET EOSCore)
    message(FATAL_ERROR "CppEOSDemo: EOSCore target not found. Build src/ first.")
endif()
if(NOT TARGET EquationOfStateV1_lib) # For includes
    message(FATAL_ERROR "CppEOSDemo: EquationOfStateV1_lib target not found for includes.")
endif()

add_executable(CppEOSDemo main.cpp)

# Link against the imported target from EOSCore
# If you used NAMESPACE EOS:: in install(EXPORT ...), it's EOS::EOSCore
target_link_libraries(CppEOSDemo PRIVATE EOSCore)