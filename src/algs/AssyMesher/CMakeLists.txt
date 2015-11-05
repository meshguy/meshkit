set(AssyMesher_srcs
  AssyMesher.cpp )

set(AssyMesher_headers
  meshkit/AssyMesher.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_library(AssyMesher STATIC
  ${AssyMesher_srcs}
  ${AssyMesher_headers})
target_link_libraries(AssyMesher
  MeshKitcore)
target_include_directories(AssyMesher
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}")
