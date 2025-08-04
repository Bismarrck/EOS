# src/cpp/math/CMakeLists.txt
add_library(EOSMath_lib STATIC
        NumericalSolvers.cpp
)

set_target_properties(EOSMath_lib PROPERTIES
        POSITION_INDEPENDENT_CODE ON
)

# Consumers of this library will need to find NumericalSolvers.h
# and will also have a dependency on Boost headers.
target_include_directories(EOSMath_lib
        PUBLIC
        # For users to #include <NumericalSolvers.h> if math/ is in include path
        # Or #include <math/NumericalSolvers.h> if src/cpp/ is in include path
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        # The install path for this header would be something like include/EOSCore/math/
        $<INSTALL_INTERFACE:include/EOSCore/math>
)

# Add installation rule for the header if it's considered public
install(FILES NumericalSolvers.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/EOSCore/math
        COMPONENT DevelopmentHeaders
)

# Add this new static library to the main export set so its dependencies
# can be tracked if necessary for installation.
install(TARGETS EOSMath_lib
        EXPORT EOSCoreTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT DevelopmentLibraries
)