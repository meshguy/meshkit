set(AssyGen_srcs
  meshscript.cpp
  AssyGen.cpp )

set(AssyGen_headers
  meshkit/AssyGen.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_library(AssyGen STATIC
  ${AssyGen_srcs}
  ${AssyGen_headers})
target_link_libraries(AssyGen
  MeshKitcore)
target_include_directories(AssyGen
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}")
