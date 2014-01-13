# /source/fem/CMakeList.txt
# Creates FEM related libraries and functions

# Define the compiler to utilize
set(CMAKE_CXX_COMPILER "mpicxx")

# Add FEMcommon library 
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/common)

# Add FEMheateq library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/heat_eq)

# Add FEMvolavg library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/volume_average)

# Build example functions
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
