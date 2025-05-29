# src/cpp/CMakeLists.txt

# For now, we will compile EquationOfStateV1.cpp as part of the test executable.
# Later, this will be:
# add_library(eos_cpp_lib STATIC
#     EquationOfStateV1.cpp
# )
# target_include_directories(eos_cpp_lib PUBLIC
#     ${CMAKE_CURRENT_SOURCE_DIR} # To find EquationOfStateV1.h from EquationOfStateV1.cpp
#     ${PROJECT_SOURCE_DIR}/include # If public headers are installed there during configure
# )

# For Phase 1, if main_phase1.cpp is in src/cpp/, we can list EquationOfStateV1.cpp
# directly in the add_executable in the root CMakeLists.txt, or create an object library
# for EquationOfStateV1.cpp here.

add_library(eos_v1_obj OBJECT
        EquationOfStateV1.cpp
)
# Ensure EquationOfStateV1.h can be found by EquationOfStateV1.cpp
target_include_directories(eos_v1_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})