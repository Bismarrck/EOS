cmake_minimum_required(VERSION 3.15)
project(EquationOfState LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Add src subdirectories
add_subdirectory(src/fortran) # Process fortran CMakeLists first
add_subdirectory(src/cpp)     # Process cpp CMakeLists

# --- Phase 1: Executable for testing ---
# This will evolve into building a library target in later phases.
add_executable(Phase1Test src/cpp/main_phase1.cpp $<TARGET_OBJECTS:eos_v1_obj>)

# Link the Fortran object files and C++ object files
target_link_libraries(Phase1Test PRIVATE
        eos_fortran_core_obj  # Fortran object files
        # eos_v1_obj # If main_phase1.cpp does not include EquationOfStateV1.cpp directly and eos_v1_obj is an OBJECT lib
        # The $<TARGET_OBJECTS:eos_v1_obj> in add_executable handles this for OBJECT libs
)

# Ensure main_phase1.cpp can find EquationOfStateV1.h
target_include_directories(Phase1Test PRIVATE ${PROJECT_SOURCE_DIR}/src/cpp)
# ------------------------------------------------------

# Optional: Add examples and tests later
# option(BUILD_EXAMPLES "Build example programs" ON)
# option(BUILD_TESTS "Build test programs" ON)

# if(BUILD_EXAMPLES)
#     add_subdirectory(examples)
# endif()

# if(BUILD_TESTS)
#     enable_testing()
#     add_subdirectory(tests)
# endif()

# Installation rules will be added later