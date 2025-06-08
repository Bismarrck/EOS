# examples/fortran_example/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(FortranEOSDemo LANGUAGES Fortran)

# Find the installed EOSCore package to get library path and C/C++ includes for the C-API
# that module_eos_fortran_interface might effectively depend on (though it's indirect).
find_package(EOSCore 1.0.0 REQUIRED)

add_executable(FortranEOSDemo main.f90)
target_sources(FortranEOSDemo PRIVATE
        main.f90
)

# Link against the imported target EOS::EOSCore
target_link_libraries(FortranEOSDemo PRIVATE eos_fortran_core_obj)
target_link_libraries(FortranEOSDemo PRIVATE EOS::EOSCore)

# Ensure the Fortran compiler can find the C headers for BIND(C) (usually implicit)
# And that it can find its own .mod files (module_eos_fortran_interface.mod)
target_include_directories(FortranEOSDemo PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR} # For .mod files generated here
)