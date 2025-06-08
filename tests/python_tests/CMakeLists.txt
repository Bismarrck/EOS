if(NOT TARGET PyEOS)
    message(WARNING "PyEOS target not found. Python tests may not be runnable via CTest without it.")
    return()
endif()

find_package(PythonInterp REQUIRED)

find_program(PYTEST_EXECUTABLE pytest)
if(PYTEST_EXECUTABLE)
    add_test(
            NAME EosPytestNumerical
            COMMAND ${PYTEST_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_pyeos_numerical.py
            --rootdir ${CMAKE_CURRENT_SOURCE_DIR} # Helps pytest find tests
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python # Ensure PyEOS is importable
    )
    # You might need to set PYTHONPATH if PyEOS isn't installed system-wide
    # set_property(TEST EosPytestNumerical PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
    set_tests_properties(EosPythonNumericalTest PROPERTIES LABELS "python_build_tree")
else()
    message(STATUS "pytest not found, skipping pytest-based Python tests.")
endif()