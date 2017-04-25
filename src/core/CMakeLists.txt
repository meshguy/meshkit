include_directories( ${MOAB_PACKAGE_INCLUDES} ${MOAB_INCLUDE_DIRS} ${CGM_INCLUDES} )

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")

set(MeshKitcore_srcs
  Error.cpp
  MeshOp.cpp
  MeshOpSet.cpp
  MeshOpSet.hpp
  MeshScheme.cpp
  MKCore.cpp
  MKGraph.cpp
  ModelEnt.cpp
  SizingFunction.cpp
  SizingFunctionVar.cpp
  VecUtil.cpp)

set(MeshKitcore_headers
  meshkit/Types.hpp
  meshkit/Error.hpp
  meshkit/GraphNode.hpp
  meshkit/MeshOp.hpp
  meshkit/MeshOpProxy.hpp
  meshkit/MeshScheme.hpp
  meshkit/MKCore.hpp
  meshkit/MKGraph.hpp
  meshkit/ModelEnt.hpp
  meshkit/NoOp.hpp
  meshkit/SizingFunction.hpp
  meshkit/SizingFunctionVar.hpp
  meshkit/VecUtil.hpp
  meshkit/Matrix.hpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/lemon
    ${CMAKE_SOURCE_DIR}/src/algs    
    ${IGEOM_INCLUDE_DIR}
    ${IGEOM_INCLUDE_DIR}/cgm)

add_library(MeshKitcore
  ${MeshKitcore_srcs}
  ${MeshKitcore_headers})

link_libraries(MeshKitcore
  MeshKitalgs
  MeshKitextern
  MeshKitutils
  iGeom
  lemon)

install(
  TARGETS   "MeshKitcore"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)


