include_directories( ${MOAB_PACKAGE_INCLUDES} ${MOAB_INCLUDE_DIRS}  ${CGM_INCLUDES} )


add_subdirectory(Qslim)
add_subdirectory(AssyGen)
add_subdirectory(AssyMesher)
add_subdirectory(CoreGen)
add_subdirectory(Sweep)
add_subdirectory(PostBL)
add_subdirectory(AdvFront)

set(extra_links)
if (ENABLE_QUADMESHER)
  add_subdirectory(QuadMesher)
  list(APPEND extra_links
    QuadMesher)
endif ()
if (ENABLE_INTASSIGN)
  add_subdirectory(IntervalAssignment)
  list(APPEND extra_links
    IntervalAssignment)
endif ()
#if (TARGET dagmc) # From MOAB
#  add_subdirectory(make_watertight)
#  list(APPEND extra_links
#    MakeWatertight)
#endif ()

set(MeshKitalgs_srcs
  register_algs.cpp
  EdgeMesher.cpp
  CurveFacetMeshReader.cpp
  SurfaceFacetMeshReader.cpp
  VertexMesher.cpp
  SCDMesh.cpp
  EBMesher.cpp
  CopyMesh.cpp
  CopyGeom.cpp
  MeshOpTemplate.cpp
  MergeMesh.cpp
  ExtrudeMesh.cpp
  CESets.cpp
  Transform.cpp)

set(MeshKitalgs_headers
  meshkit/EdgeMesher.hpp
  meshkit/CurveFacetMeshReader.hpp
  meshkit/SurfaceFacetMeshReader.hpp
  meshkit/VertexMesher.hpp
  meshkit/SCDMesh.hpp
  meshkit/EBMesher.hpp
  meshkit/CopyMesh.hpp
  meshkit/CopyGeom.hpp
  meshkit/MeshOpTemplate.hpp
  meshkit/MergeMesh.hpp
  meshkit/ExtrudeMesh.hpp
  meshkit/CESets.hpp
  meshkit/TransformBase.hpp
  meshkit/Transform.hpp)

if (MOAB_HAS_FBIGEOM)
  list(APPEND MeshKitalgs_srcs
    MBGeomOp.cpp
    MBSplitOp.cpp
    MBVolOp.cpp)
  list(APPEND MeshKitalgs_headers
    meshkit/MBGeomOp.hpp
    meshkit/MBSplitOp.hpp
    meshkit/MBVolOp.hpp)
endif ()

if (HAVE_LPSOLVER)
  add_definitions(-DHAVE_LPSOLVER)
endif ()

if (ENABLE_PARALLEL AND MOAB_HAS_PARALLEL AND CGM_HAS_PARALLEL)
  list(APPEND MeshKitalgs_srcs
    ParallelMesher.cpp
    ParExchangeMesh.cpp
    ParSendPostSurfMesh.cpp
    ParRecvSurfMesh.cpp)
  list(APPEND MeshKitalgs_headers
    meshkit/ParallelMesher.hpp
    meshkit/ParExchangeMesh.hpp
    meshkit/ParSendPostSurfMesh.hpp
    meshkit/ParRecvSurfMesh.hpp)
endif ()

add_library(MeshKitalgs
  ${MeshKitalgs_srcs}
  ${MeshKitalgs_headers})

include(GenerateExportHeader)
generate_export_header(MeshKitalgs)

link_libraries(MeshKitalgs
  MeshKitcore
  Qslim
  AssyGen
  AssyMesher
  CoreGen
  Sweep
  PostBL
  AdvFront
  ${extra_links})

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/lemon 
    ${CMAKE_SOURCE_DIR}/src/algs/Sweep    
    ${CMAKE_SOURCE_DIR}/src/algs/Qslim  
    ${CMAKE_SOURCE_DIR}/src/algs/AdvFront    
    ${CMAKE_SOURCE_DIR}/src/algs/AssyGen       
    ${CMAKE_SOURCE_DIR}/src/algs/AssyMesher
    ${CMAKE_SOURCE_DIR}/src/algs/CoreGen
    ${CMAKE_SOURCE_DIR}/src/algs/PostBL    
    )

install(
  TARGETS   "MeshKitalgs"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
