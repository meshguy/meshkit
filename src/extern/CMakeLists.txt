include_directories( ${MOAB_PACKAGE_INCLUDES} ${MOAB_INCLUDE_DIRS}  ${CGM_INCLUDES} )

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

add_library(MeshKitextern 
  STATIC register_extern.cpp)

link_libraries(MeshKitutils
    ${extra_links})

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/lemon    
    )

install(
  TARGETS   MeshKitextern
  EXPORT    MeshKit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
