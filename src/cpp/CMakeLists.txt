# src/cpp/CMakeLists.txt

add_subdirectory(math)

# Create a static library from EquationOfStateV1.cpp/.h
add_library(EquationOfStateV1_lib STATIC
        EquationOfStateV1.cpp
        eos_c_api.cpp
        materials/ComplicatedLegacyEOS.cpp
        materials/PolyEOS.cpp
        materials/AnalyticEOS.cpp
        utils/string_utils.cpp
        utils/checksum_utils.cpp
        utils/parser_utils.cpp
        utils/hdf5_utils.cpp
)

set_target_properties(EquationOfStateV1_lib PROPERTIES
        POSITION_INDEPENDENT_CODE ON
)

# Ensure EquationOfStateV1.h can be found by EquationOfStateV1.cpp
# and by consumers of this library (like tests)
target_include_directories(EquationOfStateV1_lib
        PUBLIC # For its own headers
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/utils>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        PRIVATE
        ${HDF5_INCLUDE_DIRS}
)

# Determine correct HDF5 library linking
# Modern HDF5 provides imported targets like HDF5::hdf5 and HDF5::hdf5_hl
if(TARGET HDF5::hdf5 AND TARGET HDF5::hdf5_hl)
    target_link_libraries(EquationOfStateV1_lib PRIVATE HDF5::hdf5 HDF5::hdf5_hl)
elseif(TARGET HDF5::HDF5) # Older imported target name
    target_link_libraries(EquationOfStateV1_lib PRIVATE HDF5::HDF5)
else() # Fallback to variables
    target_link_libraries(EquationOfStateV1_lib PRIVATE ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif()

target_link_libraries(EquationOfStateV1_lib PRIVATE EOSMath_lib)

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