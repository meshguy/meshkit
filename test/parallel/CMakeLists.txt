# Build and run all tests

include_directories(
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/algs
  ${CMAKE_SOURCE_DIR}/src/utils
  ${CMAKE_SOURCE_DIR}/src/lemon
  ${CMAKE_SOURCE_DIR}/src/algs/CoreGen
  ${MOAB_PACKAGE_INCLUDES}
  ${MOAB_INCLUDE_DIRS}
  ${CGM_INCLUDES}
)
include_directories(${CMAKE_SOURCE_DIR}/test)

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

set( TESTS test_coregen.cpp
    )

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" base ${fname} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "-DTEST -DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}" )
  target_link_libraries(${base}
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
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

