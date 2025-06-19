# tests/cpp_tests/CMakeLists.txt

# Create a single executable for all C++ doctests
add_executable(eos_cpp_tests
        test_main.cpp
        test_core_logic.cpp
        test_pack_unpack.cpp
        test_numerical_validation.cpp
)

# Link against the EquationOfStateV1 library and Fortran objects
target_link_libraries(eos_cpp_tests PRIVATE
        EquationOfStateV1_lib      # The C++ static library
        eos_fortran_core_obj   # The Fortran object library (defined in src/fortran/CMakeLists.txt)
)

target_include_directories(eos_cpp_tests PRIVATE
        ${HDF5_INCLUDE_DIRS}
)

# Ensure tests can find EquationOfStateV1.h (via EquationOfStateV1_lib's PUBLIC include)
# and doctest.h (via root CMakeLists.txt's include_directories)

# Add test to CTest
add_test(NAME EosCppTests COMMAND eos_cpp_tests)
set_tests_properties(EosCppTests PROPERTIES LABELS "cpp")
