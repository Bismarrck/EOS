# python/CMakeLists.txt
if(NOT TARGET EquationOfStateV1_lib OR NOT TARGET eos_fortran_core_obj)
    message(FATAL_ERROR "EquationOfStateV1_lib or eos_fortran_core_obj target not found. Ensure src/ is built first.")
endif()

# Find pybind11 using its own tools (assuming third_party/pybind11/tools is present)
# Adjust path to pybind11 if it's elsewhere relative to this CMakeLists.txt
# The THIRDPARTY_DIR is defined in the root CMakeLists.txt
include(${THIRDPARTY_DIR}/pybind11/tools/pybind11Tools.cmake)
# Alternatively, if you have pybind11 installed or use find_package:
# find_package(pybind11 REQUIRED)

# Define the Python extension module
# The output file will be <name>.so or <name>.pyd
pybind11_add_module(PyEOS SHARED
        pybind_eos.cpp
)

if(TARGET EOSCore)
    target_link_libraries(PyEOS PRIVATE EOSCore) # Link against the shared EOSCore library
else()
    # Fallback if EOSCore target isn't defined (e.g. building python module standalone)
    # This was the old way, linking components directly.
    message(STATUS "EOSCore target not found, linking PyEOS to components.")
    target_link_libraries(PyEOS PRIVATE
            EquationOfStateV1_lib
            $<TARGET_OBJECTS:eos_fortran_core_obj>
    )
endif()

# Optional: Install the Python module
# This destination should allow Python to find it, e.g. site-packages or a path added to PYTHONPATH
# pybind11_install_module(PyEOS DESTINATION ${CMAKE_INSTALL_PREFIX}/python/PyEOS COMPONENT python)
# A common pattern for local install:
install(TARGETS PyEOS
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/python_packages/PyEOS # Or some other relative path
        COMPONENT python_bindings
)