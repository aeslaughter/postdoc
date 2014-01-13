# /source/fem/examples/CMakeList.txt
# Builds the examples for the FEM library

# Create a list of the required packages
set(EXAMPLE_LIBS	
	${LIBMESH_LIBRARIES} 
	common FEMcommon
	${PETSC_LIBRARIES} 
	vtkIO vtkCommon vtkFiltering
	tbb tbbmalloc gfortran
	boost_program_options)

message(${EXAMPLE_LIBS})
include_directories(${LIBMESH_INCLUDES}) 

# Example 1
add_executable(example1 example1.cpp)
target_link_libraries(example1 FEMheateq ${EXAMPLE_LIBS})
	
# Example 2
add_executable(example2 example2.cpp)
target_link_libraries(example2 FEMheateq ${EXAMPLE_LIBS})
			
# Example 3
add_executable(example3 example3.cpp)
target_link_libraries(example3 ${EXAMPLE_LIBS})

add_executable(example3b example3b.cpp)
target_link_libraries(example3b ${EXAMPLE_LIBS})

add_executable(example3c example3c.cpp)
target_link_libraries(example3c ${EXAMPLE_LIBS})
			
# Example 4
add_executable(example4 example4.cpp)
target_link_libraries(example4 ${EXAMPLE_LIBS})

# Example 5
add_executable(example5 example5.cpp)
target_link_libraries(example5 ${EXAMPLE_LIBS})

# Example 6
add_executable(example6 example6.cpp)
target_link_libraries(example6 FEMvolavg ${EXAMPLE_LIBS})

# Test of MyMesh class
add_executable(test_mymesh test_mymesh.cpp)
target_link_libraries(test_mymesh ${EXAMPLE_LIBS})

# Test function for VolumeAverageEqData class
add_executable(test_eqdata test_eqdata.cpp)
target_link_libraries(test_eqdata ${EXAMPLE_LIBS})

# Test function for VolumeAverageEqData class
add_executable(test_thermo test_thermo.cpp)
target_link_libraries(test_thermo ${EXAMPLE_LIBS})
