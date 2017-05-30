# Build and run all tests
add_definitions("-DDATADIR=${CMAKE_SOURCE_DIR}/data")

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

add_subdirectory(core)
add_subdirectory(algs)
#add_subdirectory(extern)
if( MESHKIT_HAVE_MPI )
  add_subdirectory(parallel)
endif()

set( TESTS test_testutil.cpp)

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

