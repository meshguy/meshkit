# This makefile is still a work in progress, it has been shown
# to work only under the simplest configure with the basic example.

# MESHKIT_DIR points to top-level install dir, below which MESHKIT's lib/ and include/ are located
include ${MESHKIT_DIR}/lib/meshkit.make

.SUFFIXES: .o .cpp 

# MESH_DIR is the top-level MESHKIT source directory, used to locate mesh files that come with MESHKIT source
MESH_DIR="../data"

EXAMPLES = example_assygen example_basic example_copygeom example_copymesh example_ebmesh example_edgemesher example_extrudemesh example_fbgeom example_graph example_intassign example_jaalquad example_mbgeomop example_mbsplitop example_mbvolop example_meshoptemplate example_mesquiteopt example_ngtetmesher example_onetooneswept example_parallelmesher example_camalpaver example_postbl example_qslimmesher example_scdmesh example_setpnt2quad example_tfimapping  example_trianglemesher example_coregen example_volice

MOSTLY_CLEAN_FILES = AfterTFISetup.eps AfterWinslow.vtk InitialMapping.vtk mesquite.gpt tfimapping_out.exo

all:: $(EXAMPLES)

$(EXAMPLES): %: %.o
	${MESHKIT_CXX} -o $@ $< ${MESHKIT_LIBS_LINK}

.cpp.o :
	${MESHKIT_CXX} ${MESHKIT_CXXFLAGS} ${MESHKIT_CPP} ${MESHKIT_CPPFLAGS} ${MESHKIT_INCLUDES} ${MESHKIT_EXTERNAL_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

clean:
	rm ${MOSTLY_CLEAN_FILES} ${EXAMPLES:%=%.o} ${EXAMPLES} -f

