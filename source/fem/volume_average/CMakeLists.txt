# /source/fem/volume_average/CMakeList.txt
# FEMvolavg library

# Define the compiler to utilize
set(CMAKE_CXX_COMPILER "mpicxx")

# List of source files
set(FEM_VOLAVG_SOURCES
	#${CMAKE_CURRENT_SOURCE_DIR}/energy_system.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/momentum_system.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/concentration_system.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/thermo_system.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/level_set_system.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/front_velocity_system.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/front_velocity_eq.cpp
	)

# List the FEM dependencies
set(FEM_LIBRARIES common FEMcommon
	${LIBMESH_LIBRARIES} ${PETSC_LIBRARIES} 
	vtkIO vtkCommon vtkFiltering
	tbb tbbmalloc gfortran)	

# Add the library
add_library(FEMvolavg ${FEM_VOLAVG_SOURCES})
target_link_libraries(FEMvolavg ${FEM_VOLAVG_LIBRARIES})
