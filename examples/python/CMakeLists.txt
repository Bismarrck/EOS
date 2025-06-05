find_package(PythonInterp REQUIRED)

if(TARGET PyEOS)
    add_custom_target(RunPythonEOSTest
            COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/main.py
            DEPENDS PyEOS # Ensure Python module is built first
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} # Or where the script expects to find data/module
            COMMENT "Running Python EOS example"
            USES_TERMINAL # So output is visible
    )
    if(EOS_BUILD_TESTS) # Add to CTest if tests are enabled
        add_test(NAME EosPythonExample COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/main.py)
        set_tests_properties(EosPythonExample PROPERTIES WORKING_DIRECTORY ${CMAKE_BINARY_DIR}) # Adjust CWD if needed
    endif()
else()
    message(WARNING "Python module target PyEOS not found. Skipping RunPythonEOSTest target.")
endif()