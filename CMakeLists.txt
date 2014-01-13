# CMakeLists.txt
#-----------------------------------------------------------------------
#
# Required Environmental Variables:
# 	CGAL_DIR = /home/slaughter/packages/CGAL-4.0/; export CGAL_DIR
#	METHOD=dbg; export METHOD
#	LIBMESH_DIR=/home/slaughter/packages/libmesh-trunk; export LIBMESH_DIR
#	LIBMESH_ARCH=x86_64-unknown-linux-gnu; export LIBMESH_ARCH
#	PETSC_ARCH=linux64-gcc; export PETSC_ARCH
#	PETSC_DIR=/home/slaughter/packages/petsc-3.1-p8; export PETSC_DIR
#	MPIHOME=/usr/lib/mpich2/lib; export MPIHOME
#	VTK_DIR=/home/slaughter/packages/VTK/build; export VTK_DIR
#   
# List of Dependencies:
#	ImageMagick: libmagick++-dev must be installed, this should be done
#		using Synaptic Package Manager. This installs the needed files
#		in /usr/include/ImageMagick directory. This will need to be 
#		modified if ImageMagick is installed elsewhere.
#
#	BOOST:
#		libboost1.46-dev(all)	
#
#
# Directory Structure:
#	bin: contains all executables
#	
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# CMAKE INITILIZATION

# Define project and minimum required Cmake version
project(slaughter)
cmake_minimum_required(VERSION 2.8.4)

# Set the location for executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Set the location for the libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# Set the debug and optimized compiler flags (from libmesh make.common)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -felide-constructors -g -ansi -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -fPIC -Woverloaded-virtual -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC")

#set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")

# Toggle the build type base on the $METHOD environmental variable
if($ENV{METHOD} STREQUAL "dbg")
	set(CMAKE_BUILD_TYPE Debug)
else()
	set(CMAKE_BUILD_TYPE Release)
endif()

# Setup a custom command for generating documentation: "make doc"
add_custom_target(doc
	doxygen Doxyfile.local
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "Generating documentation with Doxygen" VERBATIM)
	
add_custom_target(web
	doxygen Doxyfile.web
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "Generating documentation for website with Doxygen" VERBATIM)
		
add_custom_target(matlab-doc
	matlab -nosplash -nodesktop -r "cd 'matlab';pub;quit;"
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "Generating MATLAB html documentation" VERBATIM)		
		
#-----------------------------------------------------------------------	
# EXTERNAL LIBRARIES

# MPICH2
include_directories($ENV{MPI_INCLUDE})
link_directories($ENV{MPIHOME})

# Find the CGAL package, the include statement creates the following 
# cmake variables that are used here:
# ${CGAL_LIBRARIES} = CGAL related libraries
# ${CGAL_3RD_PARTY_LIBRARIES} = CGAL dependencies libraries
set(CGAL_DIR $ENV{CGAL_DIR})
find_package(CGAL COMPONENTS ImageIO)
include(${CGAL_USE_FILE}) 

# Find the VTK package, the include statement adds the necessary
# library directories. Hence, only the library name,  e.g. vtkIO, is 
# requied in the target_link_libraries commands
set(VTK_DIR $ENV{VTK_DIR})
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Find ITK.
set(ITK_DIR $ENV{ITK_DIR})
find_package(ITK REQUIRED)
IF(ITK_FOUND)
  include(${ITK_USE_FILE})
  include_directories(${ITK_DIR}/Modules/Bridge/VtkGlue/include)
ENDIF(ITK_FOUND)

# Defines a list of ImageMagick, Magick++ libraries and includes the 
# necessary headers, it assumes that it was installed with Synaptic
include_directories(/usr/include/ImageMagick)
set(MAGICK_LIBRARIES Magick++ MagickCore)

# Uses the custom libMeshConfig.cmake file for incorporationg libMesh, 
# this file creates the LIBMESH_LIBRARIES and LIBMESH_INCLUDES Cmake
# variables
find_package(libMesh HINTS ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${LIBMESH_INCLUDES}) 

# PETSc Libraries
# (see www.mcs.anl.gov/petsc/documentation/faq.html#cmake) 
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)
find_package(PETSc COMPONENTS CXX REQUIRED)
add_definitions(${PETSC_DEFINITIONS})
include_directories(${PETSC_INCLUDES})

# Define BOOST related variables, this is done directly here rather than
# with the find_package command because boost was installed in a 
# non-standard location (this needs some work)
#set(BOOST_DIR /home/slaughter/packages/boost_1_49_0)
#set(BOOST_LIB_DIRS 
#    ${BOOST_DIR}/bin.v2/libs/program_options/build/gcc-4.6.1/release/threading-multi)
#include_directories(${BOOST_DIR})
#link_directories(${BOOST_LIB_DIRS})

#-----------------------------------------------------------------------
# MY LIBRARIES AND EXECUTALBES

# Adds the include directory associated with my source code
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Common library:
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/source/common)

# FEM libraries:
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/source/fem)

# Vol2mesh library:
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/source/vol2mesh)

#-----------------------------------------------------------------------
# MY TEST FUNCTIONS

add_executable(test_libmesh ${CMAKE_CURRENT_SOURCE_DIR}/source/test_libmesh.cpp)
target_link_libraries(test_libmesh
	${LIBMESH_LIBRARIES} 
	${PETSC_LIBRARIES} 
	vtkIO vtkCommon vtkFiltering
	tbb tbbmalloc gfortran)
	
# Generate the executable to test the file_parts class
#include(${CMAKE_SOURCE_DIR}/source/test/file_parts/CMakeLists.txt)
add_executable(test ${CMAKE_CURRENT_SOURCE_DIR}/source/test.cpp)
