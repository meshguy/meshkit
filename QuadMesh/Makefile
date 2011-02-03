EXT_SOFT =  /home/csverma/Software
MESQUITE_DIR=$(EXT_SOFT)/Mesh/ITAPS/Mesquite
VERDICT_DIR =$(EXT_SOFT)/Mesh/Verdict

MAIN = matching.o Tri2Quad.o
MAIN = meshclean.o Tri2Quad.o
OBJS = edmond_matching.o mesh_io.o Mesh.o DualGraph.o BinaryTree.o off_io.o laplace.o circumcenter.o QuadCleanUp.o geomPredicates.o OneDefectRemeshing.o 

OBJS = $(MAIN)  AdfrontClean.o SwapQuadEdge.o  QuadCleanUp.o OneDefectRemeshing.o DijkstraShortestPath.o  MeshRefine2D.o laplace.o

MeshObjs =  Mesh.o mesh_io.o off_io.o circumcenter.o geomPredicates.o SwapTriEdge.o tfiblend.o DualGraph.o BinaryTree.o Tri2Quad.o  QuadPatches.o MeshOpt.o Doublet.o Bridges.o Singlet.o SwapQuadEdge.o FaceClose.o QuadQuality.o 

CXX= g++ -pg
#CXX=icpc -Wall

CPPFLAGS  = -DUSE_MESQUITE -I$(MESQUITE_DIR)/include  -DUSE_VERDICT -I$(VERDICT_DIR)/include -DUSE_HASHMAP

all: meshlib fire

LIBS = libjaal.a -L$(MESQUITE_DIR)/lib -lmesquite -L$(VERDICT_DIR)/lib -lverdict110

meshlib:$(MeshObjs)
	ar rc libjaal.a $(MeshObjs)
	ranlib libjaal.a 

fire:$(OBJS) 
	$(CXX) -o fire $(OBJS) $(LIBS)

.o:.cc
	$(CXX) $(CPPFLAGS) $<

clean:
	\rm -rf *.o fire libjaal.a
