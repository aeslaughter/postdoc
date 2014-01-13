# /source/fem/common/CMakeList.txt
# FEMcommon library

# Define the compiler to utilize
set(CMAKE_CXX_COMPILER "mpicxx")

# List of source files
set(FEM_COMMON_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/my_analytic_function.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/my_mesh.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/system_base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/explicit_system_base.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/implicit_system_base.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/init_func_base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/boundary_base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/my_vtkio.cpp
)

# List the FEM dependencies
set(FEM_COMMON_LIBRARIES common
	${LIBMESH_LIBRARIES} 
	${PETSC_LIBRARIES} 
	vtkIO vtkCommon vtkFiltering
	tbb tbbmalloc gfortran)	

# Add the library
add_library(FEMcommon ${FEM_COMMON_SOURCES})
target_link_libraries(FEMcommon ${FEM_COMMON_LIBRARIES})
