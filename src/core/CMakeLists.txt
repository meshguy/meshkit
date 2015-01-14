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

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_library(MeshKitcore STATIC
  ${MeshKitcore_srcs}
  ${MeshKitcore_headers})
target_link_libraries(MeshKitcore
  MeshKitalgs
  MeshKitextern
  MeshKitutils
  iGeom
  emon)
target_include_directories(MeshKitcore
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}")
