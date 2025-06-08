# examples/example/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(CppEOSDemo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the installed EOSCore package
# CMAKE_PREFIX_PATH might need to be set by the CI script or user
# when building this example externally.
find_package(EOSCore 1.0.0 REQUIRED) # Request version 1.0.0 or newer

add_executable(CppEOSDemo main.cpp)

# Link against the imported target from EOSCore
# If you used NAMESPACE EOS:: in install(EXPORT ...), it's EOS::EOSCore
target_link_libraries(CppEOSDemo PRIVATE EOS::EOSCore)
# If no namespace, it's just EOSCore:
# target_link_libraries(CppEOSDemo PRIVATE EOSCore)

# The include directories should be handled by find_package and target_link_libraries
# if EOSCoreTargets.cmake correctly lists INTERFACE_INCLUDE_DIRECTORIES.