
#project(tools)
option ( MESHKIT_BUILD_ASSYGEN    "Build the AssyGen tool?"           ON )
option ( MESHKIT_BUILD_COREGEN    "Build the CoreGen tool?"           ON )
option ( MESHKIT_BUILD_POSTBL     "Build the PostBL tool?"           ON )
option ( MESHKIT_BUILD_MKMESH     "Build the MKmesh tool?"           ON )

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")
add_definitions("-DDATADIR=${meshkit_SOURCE_DIR}/data")

include_directories(
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/algs
  ${CMAKE_SOURCE_DIR}/src/utils
  ${CMAKE_SOURCE_DIR}/src/lemon
  ${MOAB_PACKAGE_INCLUDES}
  ${MOAB_INCLUDE_DIRS}
  ${CGM_INCLUDES}
)
include_directories(${CMAKE_SOURCE_DIR}/test)

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

link_libraries(
  MeshKitcore
  MeshKitAlgs
  CoreGen
  AssyGen
  AssyMesher
  Sweep
  PostBL
  Qslim
  AdvFront
  MeshKitextern
  MeshKitutils
  lemon
  iGeom
  MOAB
  iMesh
  iRel
)

# AssyGen
if ( MESHKIT_BUILD_ASSYGEN )
    
  add_executable(assygen
    assygen.cpp)

    include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/src/algs/CoreGen
        ${CMAKE_SOURCE_DIR}/src/algs/AssyGen
        ${CMAKE_SOURCE_DIR}/src/algs/PostBL    
        )

    set( MESHKIT_HAVE_ASSYGEN 1 )
    install(
      TARGETS assygen
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )
endif ( )
 
# AssyGen
if ( MESHKIT_BUILD_COREGEN )
    add_executable(coregen
      coregen.cpp)

    set( MESHKIT_HAVE_COREGEN 1 )
    install(
      TARGETS coregen
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )

endif ( MESHKIT_BUILD_COREGEN )

# PostBL
if ( MESHKIT_BUILD_POSTBL )
    add_executable(postbl
      postbl.cpp)

    set( MESHKIT_HAVE_POSTBL 1 )
    install(
      TARGETS postbl
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )

endif ( MESHKIT_BUILD_POSTBL )

# MKmesh
if ( MESHKIT_BUILD_MKMESH )
    add_executable(MKmesh
      MKmesh.cpp)

    set( MESHKIT_HAVE_MKMESH 1 )
    install(
      TARGETS MKmesh
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )

endif ( MESHKIT_BUILD_MKMESH )

