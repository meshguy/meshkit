# This makefile is still a work in progress, it has been shown
# to work only under the simplest configure with the basic example.

# MESHKIT_DIR points to top-level install dir, below which MESHKIT's lib/ and include/ are located
include ${MESHKIT_DIR}/lib/meshkit.make

.SUFFIXES: .o .cpp .F90

# MESH_DIR is the top-level MESHKIT source directory, used to locate mesh files that come with MESHKIT source
MESH_DIR="../data"

EXAMPLES = example_assygen example_basic example_copygeom example_copymesh example_ebmesh example_edgemesher example_extrudemesh example_fbgeom example_graph example_intassign example_jaalquad example_mb2 example_mbgeom example_MBGeomOp example_mbgeomvar example_MBSplitOp example_MBVolOp example_meshoptemplate example_mesquite example_ngtetmesher example_onetooneswept example_parallelmesher example_paver example_postbl example_qslim example_scdmesh example_SetPnt2Quad example_tetmesher example_tfimapping example_triadvance example_TriangleMesher example_volIce

MOSTLY_CLEAN_FILES = AfterTFISetup.eps AfterWinslow.vtk InitialMapping.vtk mesquite.gpt tfimapping_out.exo

all:: $(EXAMPLES)

$(EXAMPLES): %: %.o
	${MESHKIT_CXX} -o $@ $< ${MESHKIT_LIBS_LINK}

.cpp.o :
	${MESHKIT_CXX} ${MESHKIT_CXXFLAGS} ${MESHKIT_CPP} ${MESHKIT_CPPFLAGS} ${MESHKIT_INCLUDES} ${MESHKIT_EXTERNAL_INCLUDES} ${MOAB_INCLUDES} ${CGM_INT_INCLUDE} ${LASSO_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

.F90.o :
	${IMESH_FC} ${IMESH_FCFLAGS} ${IMESH_INCLUDES} ${IMESH_FCDEFS} -DMESH_DIR=\"${MESH_DIR}\" -c $<

clean:
	rm ${MOSTLY_CLEAN_FILES} ${EXAMPLES:%=%.o} ${EXAMPLES} -f

