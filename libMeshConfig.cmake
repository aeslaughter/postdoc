# libMeshConfig.cmake
# Establishes: 
#	${LIBMESH_INCLUDES}
#	${LIBMESH_LIBRARIES}
#
# Environmental varaibles:
#	METHOD=opt; export METHOD
#	LIBMESH_DIR=/home/slaughter/packages/libmesh-0.7.3.1/libmesh; export LIBMESH_DIR
#	LIBMESH_ARCH=x86_64-unknown-linux-gnu; export LIBMESH_ARCH

# Re-define libMesh related paths and triggers from environmental variables
set(LIBMESH_DIR $ENV{LIBMESH_DIR})
set(LIBMESH_ARCH_NO_METHOD $ENV{LIBMESH_ARCH})
set(LIBMESH_ARCH ${LIBMESH_ARCH_NO_METHOD}_$ENV{METHOD})

# Define variable containing libMesh library
set(LIBMESH_LIBRARIES 
	${LIBMESH_DIR}/lib/${LIBMESH_ARCH}/libmesh.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libexodusii.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libfparser.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libGK.so	
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libgmv.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libgzstream.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libHilbert.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/liblaspack.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libmetis.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libnemesis.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libnetcdf.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libparmetis.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libsfcurves.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libtetgen.so
	${LIBMESH_DIR}/contrib/lib/${LIBMESH_ARCH}/libtriangle.so	
	${LIBMESH_DIR}/contrib/tecplot/lib/${LIBMESH_ARCH_NO_METHOD}/tecio.a)

# Define libMesh include directories 
set(LIBMESH_INCLUDES
	${LIBMESH_DIR}/include
    ${LIBMESH_DIR}/include/base
    ${LIBMESH_DIR}/include/enums 
    ${LIBMESH_DIR}/include/error_estimation 
    ${LIBMESH_DIR}/include/fe
    ${LIBMESH_DIR}/include/geom 
    ${LIBMESH_DIR}/include/mesh 
    ${LIBMESH_DIR}/include/numerics 
    ${LIBMESH_DIR}/include/parallel
    ${LIBMESH_DIR}/include/partitioning
    ${LIBMESH_DIR}/include/physics
    ${LIBMESH_DIR}/include/quadrature
    ${LIBMESH_DIR}/include/reduced_basis
    ${LIBMESH_DIR}/include/solvers 
    ${LIBMESH_DIR}/include/systems
    ${LIBMESH_DIR}/include/utils 
	${LIBMESH_DIR}/contrib/exodusii/Lib/include
	${LIBMESH_DIR}/contrib/fparser
	${LIBMESH_DIR}/contrib/gmv
	${LIBMESH_DIR}/contrib/gzstream
	${LIBMESH_DIR}/contrib/libHilbert/include
	${LIBMESH_DIR}/contrib/laspack
	${LIBMESH_DIR}/contrib/metis/GKlib
	${LIBMESH_DIR}/contrib/metis/Lib
	${LIBMESH_DIR}/contrib/nemesis/Lib	
	${LIBMESH_DIR}/contrib/netcdf/Lib
	${LIBMESH_DIR}/contrib/parmetis/Lib
	${LIBMESH_DIR}/contrib/sfcurves
	${LIBMESH_DIR}/contrib/tetgen
	${LIBMESH_DIR}/contrib/triangle
	${LIBMESH_DIR}/contrib/tecplot/include)
