# /source/common/CMakeLists.txt
#
# Build my common library from its various source components

#-----------------------------------------------------------------------
# BUILD THE LIBRARY
#-----------------------------------------------------------------------

# Define the compiler to utilize
set(CMAKE_CXX_COMPILER "mpicxx")

# Set a variable containing a list of the source files
set(COMMON_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/file_parts/file_parts.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/user_options/user_options.cpp)
	
# Create the library	
add_library(common ${COMMON_SOURCES})


#-----------------------------------------------------------------------
# BUILD THE TEST FUNCTIONS
#-----------------------------------------------------------------------

# Create the executable
add_executable(test_file_parts 
	${CMAKE_CURRENT_SOURCE_DIR}/file_parts/test_file_parts.cpp)

# Link the executable to the necessary library
target_link_libraries(test_file_parts common)

# Create the executable
add_executable(test_user_options
	${CMAKE_CURRENT_SOURCE_DIR}/user_options/test_user_options.cpp)

# Link the executable to the necessary library
target_link_libraries(test_user_options 
	common
	boost_program_options)

 
