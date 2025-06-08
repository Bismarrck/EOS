# src/cpp/CMakeLists.txt

# Create a static library from EquationOfStateV1.cpp/.h
add_library(EquationOfStateV1_lib STATIC
        EquationOfStateV1.cpp
        eos_c_api.cpp
        utils/string_utils.cpp
)

set_target_properties(EquationOfStateV1_lib PROPERTIES
        POSITION_INDEPENDENT_CODE ON
)

# Ensure EquationOfStateV1.h can be found by EquationOfStateV1.cpp
# and by consumers of this library (like tests)
target_include_directories(EquationOfStateV1_lib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR} # For EquationOfStateV1.h
        ${CMAKE_CURRENT_SOURCE_DIR}/utils
)