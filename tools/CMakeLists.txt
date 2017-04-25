project(tools)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")
add_definitions("-DDATADIR=${meshkit_SOURCE_DIR}/data")
add_definitions(-DMOAB)

include_directories(${CGM_INCLUDES})
include_directories(${MOAB_INCLUDES})
include_directories(${CMAKE_SOURCE_DIR}/test)

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_executable(postbl
  postbl.cpp)

link_libraries(postbl
  MeshKitcore)

add_executable(assygen
  assygen.cpp)

link_libraries(assygen
  MeshKitalgs
  MeshKitcore
  iGeom
  iMesh
  iRel)

add_executable(coregen
  coregen.cpp)

link_libraries(coregen
  MeshKitalgs
  MeshKitcore
  iGeom
  iMesh
  MOAB
  iRel)

install(TARGETS coregen assygen
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
