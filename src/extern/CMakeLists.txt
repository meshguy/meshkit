include_directories( ${MOAB_INCLUDES} ${CGM_INCLUDES} )
set(extra_links)
if (ENABLE_CAMAL)
  add_subdirectory(CAMAL) # TODO: write CMake
  add_definitions(-DHAVE_CAMAL)
  list(APPEND extra_links
    MeshKitCAMAL)
endif ()
if (NETGEN_FOUND)
  add_subdirectory(netgen) # TODO: write CMake
  add_definitions(-DHAVE_NETGEN)
  list(APPEND extra_links
    netgen)
endif ()
if (TRIANGLE_FOUND)
  add_subdirectory(triangle) # TODO: write CMake
  add_definitions(-DHAVE_TRIANGLE)
  list(APPEND extra_links
    MeshKittriangle)
endif ()
if (MESQUITE_FOUND)
  add_subdirectory(mesquite) # TODO: write CMake
  add_definitions(-DHAVE_MESQUITE)
  list(APPEND extra_links
    MeshKitmesquite)
endif ()

add_library(MeshKitextern STATIC
  register_extern.cpp)
target_link_libraries(MeshKitextern
  MeshKitutils
  ${extra_links})
