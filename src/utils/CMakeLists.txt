set(MeshKitutils_srcs
  iMesh.cpp
  clock.cpp
  cylinder.cpp
  mstream.cpp
  parser.cpp
  pincell.cpp)

set(MeshKitutils_headers
  meshkit/iBase.hpp
  meshkit/iGeom.hpp
  meshkit/iMesh.hpp
  meshkit/iRel.hpp
  meshkit/LocalSet.hpp
  meshkit/LocalTag.hpp
  meshkit/RegisterMeshOp.hpp
  meshkit/SimpleArray.hpp
  meshkit/parser.hpp
  meshkit/pincell.hpp
  meshkit/clock.hpp
  meshkit/cylinder.hpp
  meshkit/vectortemplate.hpp
  meshkit/matrixtemplate.hpp
  meshkit/mstream.hpp
  meshkit/ReadPolyLine.hpp)

if (MOAB_HAS_FBIGEOM)
  list(APPEND MeshKitutils_headers
    meshkit/FBiGeom.hpp)
endif ()

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")
if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_library(MeshKitutils 
  ${MeshKitutils_srcs}
  ${MeshKitutils_headers})

link_libraries(MeshKitutils
  MeshKitcore)

include_directories(MeshKitutils
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CGM_INCLUDE_DIRS}
    ${MOAB_INCLUDE_DIRS}
    ${Lasso_INCLUDE_DIRS})

install(
  TARGETS   "MeshKitutils"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
