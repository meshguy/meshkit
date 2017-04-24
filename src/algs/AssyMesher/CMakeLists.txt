set(AssyMesher_srcs
  AssyMesher.cpp )

set(AssyMesher_headers
  meshkit/AssyMesher.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_library(AssyMesher
  ${AssyMesher_srcs}
  ${AssyMesher_headers})

link_libraries(AssyMesher
  MeshKitcore)

include_directories(
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

install(
  TARGETS   "AssyMesher"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
