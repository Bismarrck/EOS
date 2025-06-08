# examples/example/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(CppEOSDemo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(CppEOSDemo main.cpp)

if(EOS_LIBRARY_IS_BUILDING_EXAMPLES)
    # This example is part of the main project build.
    # Link directly to the targets defined in the parent project.
    # Ensure these targets (EOSCore, EquationOfStateV1_lib) are defined before this point.
    if(NOT TARGET EOSCore)
        message(FATAL_ERROR "CppEOSDemo: EOSCore target not found. Build src/ first.")
    endif()
    if(NOT TARGET EquationOfStateV1_lib) # For includes
        message(FATAL_ERROR "CppEOSDemo: EquationOfStateV1_lib target not found for includes.")
    endif()
    # Link against the imported target from EOSCore
    # If you used NAMESPACE EOS:: in install(EXPORT ...), it's EOS::EOSCore
    target_link_libraries(CppEOSDemo PRIVATE EOSCore)
    target_compile_definitions(CppEOSDemo PRIVATE BUILT_WITHIN_MAIN_PROJECT)
    target_include_directories(CppEOSDemo PRIVATE
            ${CMAKE_SOURCE_DIR}/src/cpp
            ${CMAKE_SOURCE_DIR}/src/cpp/utils
    )

else()
    # --- Building as a standalone example (e.g., CI testing installed library) ---
    message(STATUS "CppEOSDemo: Building against installed EOSCore (using find_package).")

    # CMAKE_PREFIX_PATH should be set by the CI script (or user) to point to the install_dir
    find_package(EOSCore 1.0.0 REQUIRED) # Finds installed library

    # Link against the imported target (e.g., EOS::EOSCore)
    target_link_libraries(CppEOSDemo PRIVATE EOS::EOSCore)

    # The include directories for <EOSCore/EquationOfStateV1.h> are provided by find_package
    # and target_link_libraries against the imported target.
    # No target_compile_definitions(CppEOSDemo PRIVATE BUILT_WITHIN_MAIN_PROJECT) here.
endif()
