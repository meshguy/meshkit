set(PostBL_srcs
  PostBL.cpp)

set(PostBL_headers
  meshkit/PostBL.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")
if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_library(PostBL STATIC
  ${PostBL_srcs}
  ${PostBL_headers})
target_link_libraries(PostBL
  MeshKitcore)
target_include_directories(PostBL
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}")
