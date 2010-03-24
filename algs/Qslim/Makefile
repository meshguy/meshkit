
CXX=/homes/sharp/3rdparty/mpich2/mpich2-1.1.1p1/gcc/bin/mpicxx

.cpp.o:
	$(CXX) $(CFLAGS) -c $<


GFX2_FLAGS =  

include /homes/iulian/lib/MOAB/lib/moab.make
include /homes/iulian/lib/MOAB/lib/iMesh-Defs.inc

# debug options
OPTFLAGS = -O0 -g $(GFX2_FLAGS)
FEATURES =

CFLAGS = $(IMESH_INCLUDES) $(MOAB_INCLUDES) $(OPTFLAGS) 
LFLAGS = $(LIBDIR) $(OPTFLAGS) 

CORE_SRCS = Qslim_test.cpp QslimOptions.cpp QslimDecimation.cpp \
              3D.cpp  quadrics.cpp \
	     Heap.cpp Mat4.cpp primitives.cpp

QSLIM_SRCS = $(CORE_SRCS)
QSLIM_OBJS = $(QSLIM_SRCS:.cpp=.o)



qslim2     : $(QSLIM_OBJS)
	$(CXX) -o $@ $(QSLIM_OBJS) -lm $(IMESH_LIBS) $(MOAB_LIBS_LINK)


clean:
	rm -rf *.o qst
