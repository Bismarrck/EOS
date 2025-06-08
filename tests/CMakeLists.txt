# tests/CMakeLists.txt
if(EOS_BUILD_TESTS) # Assuming EOS_BUILD_TESTS is defined in parent scope
    add_subdirectory(cpp_tests)
    add_subdirectory(fortran_tests)

    if(TARGET PyEOS AND EOS_BUILD_PYTHON_BINDINGS) # If Python bindings are built
        add_subdirectory(python_tests)
    endif()
endif()
