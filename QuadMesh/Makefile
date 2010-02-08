EXT_SOFT =  /u/c/s/csverma/Software
#EXT_SOFT =  /home/csverma/Software
MOAB_DIR   = $(EXT_SOFT)/Mesh/ITAPS/MOAB
#CGM_DIR    = $(EXT_SOFT)/Mesh/ITAPS/CGM
#HDF5_DIR   = $(EXT_SOFT)/Database/HDF5
#CUBIT_DIR  = $(EXT_SOFT)/Mesh/Cubit10.2
#CUBIT_LIBS   = -L$(CUBIT_DIR)/bin -lcubiti19 -lgtcAttrib
MESQUITE_DIR=$(EXT_SOFT)/Mesh/ITAPS/Mesquite


CPPFLAGS = -O3 -DUSE_MOAB -DUSE_MESQUITE -DUSE_HASHMAP -I$(EXT_SOFT)/Boost/include -I$(MESQUITE_DIR)/include -I$(MOAB_DIR)/include -I$(HDF5_DIR)/include -I$(NETCDF_DIR)/include 

#LIBS     = -L$(MOAB_DIR)/lib -liMesh -lMOAB -L$(NETCDF_DIR)/lib -lnetcdf_c++ -lnetcdf -L$(HDF5_DIR)/lib -lhdf5 $(CUBIT_LIBS)  -L$(CGM_DIR)/lib -lcgm -L$(MESQUITE_DIR)/lib -lmesquite

LIBS     = -L$(MOAB_DIR)/lib -liMesh -lMOAB  -L$(MESQUITE_DIR)/lib  -lmesquite

OBJS = matching.o tri_io.o Mesh.o DualGraph.o BinaryTree.o Tri2Quad.o off_io.o laplace.o QuadCleanUp.o 
#OBJS =  meshOpt.o

fire:$(OBJS)
	g++ -o fire $(OBJS) $(LIBS)
.o:.cc
	g++ $(CPPFLAGS) $<
clean:
	\rm -rf *.o

