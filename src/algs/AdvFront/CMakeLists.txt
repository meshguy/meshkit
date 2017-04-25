include_directories( ${MOAB_PACKAGE_INCLUDES} ${MOAB_INCLUDE_DIRS}  ${CGM_INCLUDES} )
set(AdvFront_srcs
  	AF2Algorithm.cpp 
	AF2AlgorithmResult.cpp 
	AF2Binding.cpp 
	AF2DfltPlaneProjMaker.cpp 
	AF2DfltRuleAppVisitor.cpp 
	AF2DfltTriangleRules.cpp 
	AF2DfltTriangleMeshOp.cpp 
	AF2Edge2D.cpp 
	AF2Edge3D.cpp 
	AF2FreeZone.cpp 
	AF2FreeZoneDefLCQualLim.cpp 
	AF2FreeZoneDefSimple.cpp 
	AF2Front.cpp 
	AF2Neighborhood.cpp 
	AF2PlaneProjection.cpp 
	AF2PntTrnsfrmLnrV.cpp 
	AF2Point3D.cpp 
	AF2PointTransformNone.cpp 
    AF2Polygon2D.cpp 
    AF2Polygon3D.cpp 
	AF2Rule.cpp 
    AF2RuleApplication.cpp 
	AF2RuleExistEdge.cpp 
	AF2RuleExistVertex.cpp 
	AF2RuleNewEdge.cpp 
	AF2RuleNewTriangle.cpp 
	AF2RuleNewVertex.cpp)

set(AdvFront_headers
  	meshkit/AF2Algorithm.hpp 
	meshkit/AF2AlgorithmResult.hpp 
	meshkit/AF2Binding.hpp 
	meshkit/AF2DfltPlaneProjMaker.hpp 
	meshkit/AF2DfltRuleAppVisitor.hpp 
	meshkit/AF2DfltTriangleRules.hpp 
	meshkit/AF2DfltTriangleMeshOp.hpp 
	meshkit/AF2Edge2D.hpp 
	meshkit/AF2Edge3D.hpp 
	meshkit/AF2FreeZone.hpp 
	meshkit/AF2FreeZoneDef.hpp 
	meshkit/AF2FreeZoneDefLCQualLim.hpp 
	meshkit/AF2FreeZoneDefSimple.hpp 
	meshkit/AF2Front.hpp 
	meshkit/AF2LocalTransform.hpp 
	meshkit/AF2LocalTransformMaker.hpp 
	meshkit/AF2Neighborhood.hpp 
	meshkit/AF2PlaneProjection.hpp 
	meshkit/AF2PntTrnsfrmLnrV.hpp 
	meshkit/AF2Point2D.hpp 
	meshkit/AF2Point3D.hpp 
	meshkit/AF2PointTransform.hpp 
	meshkit/AF2PointTransformNone.hpp 
    meshkit/AF2Polygon2D.hpp 
    meshkit/AF2Polygon3D.hpp 
	meshkit/AF2Rule.hpp 
    meshkit/AF2RuleApplication.hpp 
    meshkit/AF2RuleAppVisitor.hpp 
	meshkit/AF2RuleExistEdge.hpp 
	meshkit/AF2RuleExistVertex.hpp 
	meshkit/AF2RuleNewEdge.hpp 
	meshkit/AF2RuleNewFace.hpp 
	meshkit/AF2RuleNewTriangle.hpp 
	meshkit/AF2RuleNewVertex.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")
if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/lemon)

add_library(AdvFront 
  ${AdvFront_srcs}
  ${AdvFront_headers})

link_libraries(AdvFront
  MeshKitcore)

include_directories(
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

install(
  TARGETS   "AdvFront"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
