set(Sweep_srcs
  Global.hpp
  OneToOneSwept.cpp
  TFIMapping.cpp
  EquipotentialSmooth.cpp
  EquipotentialSmooth.hpp
  GreenCoordinates3D.cpp
  GreenCoordinates3D.hpp
  IsoLaplace.cpp
  IsoLaplace.hpp
  Dijkstra.hpp
  Dijkstra.cpp)
  #HarmonicMap.cpp
  #HarmonicMap.hpp
  #MatrixOp.hpp
  #MatrixOp.cpp
  #SmoothAngleBased.hpp
  #SmoothAngleBased.cpp)

set(Sweep_headers
  meshkit/OneToOneSwept.hpp
  Global.hpp
  meshkit/TFIMapping.hpp)

if (ARMADILLO_FOUND)
  add_definitions(-DHAVE_ARMADILLO)
endif ()

if (OPENBLAS_FOUND)
  add_definitions(-DHAVE_OPENBLAS)
endif ()

if (ARMADILLO_FOUND AND OPENBLAS_FOUND)
  list(APPEND Sweep_srcs
    TriharmonicRBF.hpp
    TriharmonicRBF.cpp
    HarmonicMapper.hpp
    HarmonicMapper.cpp
    SurfProHarmonicMap.hpp
    SurfProHarmonicMap.cpp
    Deform2D.hpp
    Deform2D.cpp)
endif ()

if (MESQUITE_FOUND)
  add_definitions(-DHAVE_MESQUITE)
  list(APPEND Sweep_srcs
    SweepWrapper.hpp
    SweepWrapper.cpp
    MeshImprove.cpp
    SmartLaplaceWrapper.hpp
    SmartLaplaceWrapper.cpp)
  list(APPEND Sweep_headers
    meshkit/MeshImprove.hpp)
endif ()

if (LPSOLVER_FOUND)
  add_definitions(-DHAVE_LPSOLVER)
  list(APPEND Sweep_srcs
    SubMapping.cpp
    LPSolveClass.hpp
    LPSolveClass.cpp)
  list(APPEND Sweep_headers
    meshkit/SubMapping.hpp)
endif ()

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")

add_library(Sweep 
  ${Sweep_srcs}
  ${Sweep_headers})

link_libraries(Sweep
  MeshKitalgs
  MeshKitcore
  iRel)

include_directories(
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/core    
    ${CMAKE_SOURCE_DIR}/src/lemon
    ${CMAKE_SOURCE_DIR}/src/algs
    )

install(
  TARGETS   "Sweep"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
