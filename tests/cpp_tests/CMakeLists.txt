# tests/cpp_tests/CMakeLists.txt
enable_testing() # Enables CTest integration

add_executable(EosCppTests
    test_main.cpp          # Contains main() for doctest
    test_core_logic.cpp
    test_pack_unpack.cpp
    # Add other test_*.cpp files here
)

# Link against the EOS library components (objects for now, library later)
# We need the C++ objects from EquationOfStateV1 and the Fortran objects
target_link_libraries(EosCppTests PRIVATE
    eos_v1_obj           # From src/cpp/CMakeLists.txt
    eos_fortran_core_obj # From src/fortran/CMakeLists.txt
)
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_libraries(EosCppTests PRIVATE stdc++fs) # If any part still needs it, but we removed it
endif()


# Include directories needed by the tests
target_include_directories(EosCppTests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}         # For doctest.h if it's here
    ${PROJECT_SOURCE_DIR}/src/cpp      # For EquationOfStateV1.h
)

# Add test to CTest
add_test(NAME EosCoreFunctionalityTest COMMAND EosCppTests)
# You can add more specific tests or command line args for doctest here
# e.g., add_test(NAME EosPackUnpackTest COMMAND EosCppTests -tc="*Pack/Unpack*")