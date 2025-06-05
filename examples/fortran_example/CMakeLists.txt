# examples/fortran_example/CMakeLists.txt
if(NOT TARGET EquationOfStateV1_lib OR NOT TARGET eos_fortran_core_obj)
    message(FATAL_ERROR "EquationOfStateV1_lib or eos_fortran_core_obj target not found. Build src/ first.")
endif()

add_executable(FortranEOSTest main.f90)

target_link_libraries(FortranEOSTest PRIVATE
        EquationOfStateV1_lib    # C++ static library (includes C-API code)
        eos_fortran_core_obj # Fortran object library (includes all Fortran modules)
)

get_target_property(EOS_MOD_OUT_DIR eos_fortran_core_obj BINARY_DIR)
if(EOS_MOD_OUT_DIR)
    # For object libraries, mods are often in a subdirectory of BINARY_DIR
    # e.g.CMakeFiles/eos_fortran_core_obj.dir/ if not using CMAKE_Fortran_MODULE_OUTDIR
    # This assumes a default layout or that CMAKE_Fortran_MODULE_DIRECTORY was set for the target.
    # For an OBJECT library, the .mod files are typically created alongside the .o files.
    # The most reliable way is that the eos_fortran_core_obj target itself should
    # declare its module directory as an INTERFACE property.
    #
    # In src/fortran/CMakeLists.txt:
    # target_include_directories(eos_fortran_core_obj INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
    # This makes its build directory available for .mod files.
    # Let's assume this is done.
else()
    # Fallback if property not found (should not happen for existing target)
    # message(WARNING "Could not get BINARY_DIR for eos_fortran_core_obj. .mod files might not be found.")
    # Try a common relative path:
    include_directories(${CMAKE_BINARY_DIR}/src/fortran)
endif()

# Add Fortran example to CTest
if(EOS_BUILD_TESTS)
    add_test(NAME EosFortranExample COMMAND FortranEOSTest)
endif()