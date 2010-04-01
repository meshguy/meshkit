EXT_SOFT =  /u/c/s/csverma/Software
#EXT_SOFT =  /home/csverma/Software
#MOAB_DIR   = $(EXT_SOFT)/Mesh/ITAPS/MOAB
#CGM_DIR    = $(EXT_SOFT)/Mesh/ITAPS/CGM
#HDF5_DIR   = $(EXT_SOFT)/Database/HDF5
#CUBIT_DIR  = $(EXT_SOFT)/Mesh/Cubit10.2
#CUBIT_LIBS   = -L$(CUBIT_DIR)/bin -lcubiti19 -lgtcAttrib
#BOOST_DIR  =  $(EXT_SOFT)/Boost
MESQUITE_DIR=$(EXT_SOFT)/Mesh/ITAPS/Mesquite

CPPFLAGS = -O3 -DUSE_MOAB -DUSE_MESQUITE -I$(MESQUITE_DIR)/include -I$(MOAB_DIR)/include -I$(HDF5_DIR)/include -I$(NETCDF_DIR)/include -I$(BOOST_DIR)/include 

LIBS     = -L$(MOAB_DIR)/lib -liMesh -lMOAB  -L$(MESQUITE_DIR)/lib  -lmesquite 
#-L$(HDF5_DIR)/lib -lhdf5 
#$(CUBIT_LIBS) -L$(CGM_DIR)/lib -lcgm  -lnetcdf_c++ -lnetcdf

OBJS = matching.o mesh_io.o Mesh.o DualGraph.o BinaryTree.o Tri2Quad.o off_io.o laplace.o QuadCleanUp.o  circumcenter.o geomPredicates.o MeshRefine2D.o EdgeFlip.o

#OBJS = edmond_matching.o mesh_io.o Mesh.o DualGraph.o BinaryTree.o Tri2Quad.o off_io.o laplace.o QuadCleanUp.o 

fire:$(OBJS)
	g++ -O3 -o fire $(OBJS) $(LIBS)

QGL_DIR  = $(EXT_SOFT)/Graphics/QGLViewer
QGL_LIBS = -L$(QGL_DIR)/lib -lQGLViewer -lGL -lGLU -lQtCore

DEFINES       = -DQT_NO_DEBUG -DQT_XML_LIB -DQT_OPENGL_LIB -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
CFLAGS        = -pipe -O2 -D_REENTRANT  -W $(DEFINES)
CXXFLAGS      = -pipe -O2 -D_REENTRANT  -W $(DEFINES)

QINC = $(CXXFLAGS) -I/usr/share/qt4/mkspecs/linux-g++ -I. -I/usr/include/qt4/QtCore -I/usr/include/qt4/QtGui -I/usr/include/qt4/QtOpenGL -I/usr/include/qt4/QtXml -I/usr/include/qt4 -I../.. -I/usr/X11R6/include -I.moc -I/home/csverma/Software/Boost/include -I$(QGL_DIR)/include 


viewer:meshViewer.o
	g++ -o viewer meshViewer.o $(QGL_LIBS)

meshViewer.o: meshViewer.cc
	g++ $(QINC) -c meshViewer.cc

.o:.cc
	g++ $(CPPFLAGS) $<

clean:
	\rm -rf *.o viewer
