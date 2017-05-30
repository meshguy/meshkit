set(PostBL_srcs
  createblelements.cpp
  pushbmesh.cpp
  computenormals.cpp
  PostBL.cpp)

set(PostBL_headers
  meshkit/PostBL.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")
if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_library(PostBL
  ${PostBL_srcs}
  ${PostBL_headers})

link_libraries(PostBL
  MeshKitcore)

include_directories(
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/lemon    
    )

install(
  TARGETS   PostBL
  EXPORT    MeshKitAlgs
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)

