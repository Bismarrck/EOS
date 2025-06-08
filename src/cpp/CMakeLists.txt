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
        # For builds within this project (e.g., EOSCore linking EquationOfStateV1_lib)
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/utils>
        $<INSTALL_INTERFACE:include> # This path is relative to CMAKE_INSTALL_PREFIX
)

# Install EquationOfStateV1_lib and add it to the EOSCoreTargets export set
# This assumes CMAKE_INSTALL_LIBDIR is defined (e.g., by include(GNUInstallDirs) in root)
install(TARGETS EquationOfStateV1_lib
        EXPORT EOSCoreTargets # Same export set as EOSCore
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT DevelopmentLibraries
        # For static libs, ARCHIVE is the primary output.
        # PUBLIC_HEADER can also be specified here if its headers are distinct and need different install rules
        # PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/EquationOfStateV1_lib COMPONENT DevelopmentHeaders
        # But since its headers are already installed via root for EOSCore, this might be redundant.
)