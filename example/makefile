# This makefile is still a work in progress, it has been shown
# to work only under the simplest configure with the basic example.

# MESHKIT_DIR points to top-level install dir, below which MESHKIT's lib/ and include/ are located
include ${MESHKIT_DIR}/lib/meshkit.make

.SUFFIXES: .o .cpp .F90

# MESH_DIR is the top-level MESHKIT source directory, used to locate mesh files that come with MESHKIT source
MESH_DIR="../data"

EXAMPLES = assygen basic copygeom copymesh ebmesh edgemesh extrudemesh fbgeom graph intassign jaalquad mb2 mbgeom MBGeomOp mbgeomvar MBSplitOp MBVolOp meshoptemplate mesquite ngtetmesher onetooneswept parallelmesher paver postbl qslim scdmesh SetPnt2Quad tetmesher tfimapping triadvance TriangleMesher volIce

all:: $(EXAMPLES)

$(EXAMPLES): %: %.o
	${MESHKIT_CXX} -o $@ $< ${MESHKIT_LIBS_LINK}

.cpp.o :
	${MESHKIT_CXX} ${MESHKIT_CXXFLAGS} ${MESHKIT_CPP} ${MESHKIT_CPPFLAGS} ${MESHKIT_INCLUDES} ${MESHKIT_EXTERNAL_INCLUDES} ${MOAB_INCLUDES} ${CGM_INT_INCLUDE} ${LASSO_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

.F90.o :
	${IMESH_FC} ${IMESH_FCFLAGS} ${IMESH_INCLUDES} ${IMESH_FCDEFS} -DMESH_DIR=\"${MESH_DIR}\" -c $<

clean:
	rm ${EXAMPLES:%=%.o} ${EXAMPLES} -f

