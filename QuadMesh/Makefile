EXT_SOFT =  /u/c/s/csverma/Software
#EXT_SOFT =  /home/csverma/Software
MOAB_DIR   = $(EXT_SOFT)/Mesh/ITAPS/MOAB
#CGM_DIR    = $(EXT_SOFT)/Mesh/ITAPS/CGM
#HDF5_DIR   = $(EXT_SOFT)/Database/HDF5
#CUBIT_DIR  = $(EXT_SOFT)/Mesh/Cubit10.2
#CUBIT_LIBS   = -L$(CUBIT_DIR)/bin -lcubiti19 -lgtcAttrib
BOOST_DIR  =  $(EXT_SOFT)/Boost
MESQUITE_DIR=$(EXT_SOFT)/Mesh/ITAPS/Mesquite


CPPFLAGS = -O3 -DEXACT -DUSE_MOAB -DUSE_MESQUITE -I$(MESQUITE_DIR)/include -I$(MOAB_DIR)/include -I$(HDF5_DIR)/include -I$(NETCDF_DIR)/include -I$(BOOST_DIR)/include 

LIBS     = -L$(MOAB_DIR)/lib -liMesh -lMOAB  -L$(MESQUITE_DIR)/lib  -lmesquite

OBJS = matching.o mesh_io.o Mesh.o DualGraph.o BinaryTree.o Tri2Quad.o off_io.o laplace.o QuadCleanUp.o  circumcenter.o  geomPredicates.o

#OBJS = edmond_matching.o mesh_io.o Mesh.o DualGraph.o BinaryTree.o Tri2Quad.o off_io.o laplace.o QuadCleanUp.o 

fire:$(OBJS)
	g++ -O3 -o fire $(OBJS) $(LIBS)
	cp -rf fire ../bin
.o:.cc
	g++ $(CPPFLAGS) $<
clean:
	\rm -rf *.o
