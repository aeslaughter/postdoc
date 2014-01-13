# /source/fem/heat_eq/CMakeList.txt
# FEMheateq library

# Define the compiler to utilize
set(CMAKE_CXX_COMPILER "mpicxx")

# List of source files
set(FEM_HEATEQ_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/heat_eq_boundary.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/heat_eq_system.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/heat_eq.cpp	
	)

# List the FEM dependencies
set(FEM_LIBRARIES common FEMcommon
	${LIBMESH_LIBRARIES} ${PETSC_LIBRARIES} 
	vtkIO vtkCommon vtkFiltering
	tbb tbbmalloc gfortran)	

# Add the library
add_library(FEMheateq ${FEM_HEATEQ_SOURCES})
target_link_libraries(FEMheateq ${FEM_HEATEQ_LIBRARIES})
