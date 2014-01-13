# /source/vol2mesh/CMakeLists.txt
#
# Create the vol2mesh program executable, it links the various libraries
# needed by the program and creates the executable in the bin directory
# Define project and minimum required Cmake version

# Use the standard c++ compilier
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-O2")

# Set a variable containing a list of the source files
set(VOL2MESH_SOURCES vol2mesh.cpp vol2mesh_mesh_criteria.cpp)

# Create a library containing vol2mesh classes and functions
add_library(vol2mesh ${VOL2MESH_SOURCES})
	
# Link the associated libraries
target_link_libraries(vol2mesh common
	${CGAL_LIBRARIES} 
	${CGAL_3RD_PARTY_LIBRARIES} 
	${ITK_LIBRARIES}	
	vtkHybrid vtkIO vtkexoIIc vtkCommon vtkFiltering vtkParallel)	
	
#Create the executable 	
add_executable(v2m v2m.cpp)
target_link_libraries(v2m common vol2mesh
	boost_program_options)	


