EXT_SOFT =  /home/csverma/Software
MOAB_DIR =  $(EXT_SOFT)/Mesh/ITAPS/MOAB
MESQUITE_DIR=$(EXT_SOFT)/Mesh/ITAPS/Mesquite
VERDICT_DIR =$(EXT_SOFT)/Mesh/Verdict

MAIN = meshclean.o 
MAIN = matching.o JaalMoabConverter.o
OBJS = edmond_matching.o mesh_io.o Mesh.o DualGraph.o BinaryTree.o off_io.o laplace.o circumcenter.o QuadCleanUp.o geomPredicates.o OneDefectRemeshing.o 

OBJS = $(MAIN)  

MeshObjs =  Mesh.o mesh_io.o off_io.o circumcenter.o geomPredicates.o SwapTriEdge.o tfiblend.o DualGraph.o BinaryTree.o Tri2Quad.o  QuadPatches.o MeshOpt.o Doublet.o Bridges.o Singlet.o SwapQuadEdge.o FaceClose.o QuadQuality.o AdfrontClean.o SwapQuadEdge.o  QuadCleanUp.o OneDefectRemeshing.o DijkstraShortestPath.o  MeshRefine2D.o laplace.o 

CXX= g++ -O4


CPPFLAGS  = -I. -DUSE_MOAB -I$(MOAB_DIR)/include -DUSE_MESQUITE -I$(MESQUITE_DIR)/include  -DUSE_VERDICT -I$(VERDICT_DIR)/include -DUSE_HASHMAP

all: meshlib fire

LIBS = -L$(MESQUITE_DIR)/lib -lmesquite -L$(VERDICT_DIR)/lib -lverdict110 -L$(MOAB_DIR)/lib -liMesh -lMOAB -lhdf5 -lz ./libjaal.a

meshlib:$(MeshObjs)
	ar rc libjaal.a $(MeshObjs)
	ranlib libjaal.a 

fire:$(OBJS) 
	$(CXX) -o fire $(OBJS) $(LIBS)

.o:.cc
	$(CXX) $(CPPFLAGS) $<

clean:
	\rm -rf *.o fire libjaal.a
