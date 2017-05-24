include_directories( ${MOAB_PACKAGE_INCLUDES} ${MOAB_INCLUDE_DIRS}  ${CGM_INCLUDES} )

#project(tools)
option ( MESHKIT_BUILD_ASSYGEN    "Build the AssyGen tool?"           ON )
option ( MESHKIT_BUILD_COREGEN    "Build the CoreGen tool?"           ON )
option ( MESHKIT_BUILD_POSTBL      "Build the PostBL tool?"           ON )
option ( MESHKIT_BUILD_POSTBL      "Build the MKmesh tool?"           ON )

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")
add_definitions("-DDATADIR=${meshkit_SOURCE_DIR}/data")

include_directories(${CGM_INCLUDES})
include_directories(${MOAB_INCLUDE_DIRSS} {MOAB_PACKAGE_INCLUDES} )
include_directories(${CMAKE_SOURCE_DIR}/test)

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

link_libraries(MeshKitcore CoreGen AssyGen AssyMesher Sweep PostBL Qslim AdvFront MeshKitAlgs MeshKitextern MeshKitutils lemon)

# AssyGen
if ( MESHKIT_BUILD_ASSYGEN )
    
    link_libraries(
      iGeom
      MOAB
      iMesh
      iRel
    )
  
  add_executable(assygen
    assygen.cpp)

    include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/src/utils
        ${CMAKE_SOURCE_DIR}/src/core
        ${CMAKE_SOURCE_DIR}/src/lemon    
        ${CMAKE_SOURCE_DIR}/src/algs/CoreGen
        ${CMAKE_SOURCE_DIR}/src/algs/AssyGen
        ${CMAKE_SOURCE_DIR}/src/algs/PostBL    
        )
#    set_source_files_properties( size.cpp measure.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
#    target_link_libraries( mbsize ${TOOLS_LINK_LIBS} )
    set( MESHKIT_HAVE_ASSYGEN 1 )
    install(
      TARGETS assygen
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )
endif ( )
    
#add_executable(postbl
#  postbl.cpp)

#link_libraries(postbl
#    MeshKitalgs
#    MeshKitcore
#    iGeom
#    iMesh
#    MOAB
#    iRel)

#add_executable(coregen
#  coregen.cpp)

#link_libraries(coregen
#  MeshKitalgs
#  MeshKitcore
#  iGeom
#  iMesh
#  MOAB
#  iRel)

#install(TARGETS coregen assygen
#        RUNTIME DESTINATION bin
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib)
