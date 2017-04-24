set(AssyGen_srcs
  meshscript.cpp
  AssyGen.cpp )

set(AssyGen_headers
  meshkit/AssyGen.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_library(AssyGen
  ${AssyGen_srcs}
  ${AssyGen_headers})

link_libraries(AssyGen
  MeshKitcore)

include_directories(
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

install(
  TARGETS   "AssyGen"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
