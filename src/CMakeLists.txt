# This CMakeLists.txt in src/ primarily serves to group sub-projects
# or define targets that span across cpp and fortran subdirectories.

# For Phase 0, we will compile Fortran sources into an object library.
# In later phases, this will become a proper static/shared library.
add_subdirectory(fortran)
# add_subdirectory(cpp) # We'll add this when cpp sources are part of a library