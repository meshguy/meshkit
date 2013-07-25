# This makefile is still a work in progress, it has been shown
# to work only under the simplest configure with the basic example.

# MESHKIT_DIR points to top-level install dir, below which MESHKIT's lib/ and include/ are located
include ${MESHKIT_DIR}/lib/meshkit.make
#include /homes/brsquared/lib/meshkit/lib/meshkit.make
#include ${MESHKIT_DIR}/lib/iMesh-Defs.inc

.SUFFIXES: .o .cpp .F90

# MESH_DIR is the top-level MESHKIT source directory, used to locate mesh files that come with MESHKIT source
MESH_DIR="../data"

basic : basic.o
	${MESHKIT_CXX} -o $@ $< ${MESHKIT_LIBS_LINK}

copygeom : copygeom.o
	${MESHKIT_CXX} -o $@ $< ${MESHKIT_LIBS_LINK}

postbl : postbl.o
	${MESHKIT_CXX} -o $@ $< ${MESHKIT_LIBS_LINK}

.cpp.o :
	${MESHKIT_CXX} ${MESHKIT_CXXFLAGS} ${MESHKIT_INCLUDES} ${MESHKIT_EXTERNAL_INCLUDES} ${MOAB_INCLUDES} ${CGM_INT_INCLUDE} ${LASSO_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

.F90.o :
	${IMESH_FC} ${IMESH_FCFLAGS} ${IMESH_INCLUDES} ${IMESH_FCDEFS} -DMESH_DIR=\"${MESH_DIR}\" -c $<

