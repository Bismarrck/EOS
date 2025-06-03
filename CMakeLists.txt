cmake_minimum_required(VERSION 3.15)
project(EquationOfState LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Add src subdirectories
add_subdirectory(src/fortran) # Process fortran CMakeLists first
add_subdirectory(src/cpp)     # Process cpp CMakeLists
