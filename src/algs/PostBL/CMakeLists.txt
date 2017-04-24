set(PostBL_srcs
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
    ${CMAKE_CURRENT_SOURCE_DIR})

install(
  TARGETS   "PostBL"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)

