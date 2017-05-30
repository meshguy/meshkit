# Build and run all tests

include_directories(
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/algs
  ${CMAKE_SOURCE_DIR}/src/algs/AssyGen
  ${CMAKE_SOURCE_DIR}/src/algs/PostBL
  ${CMAKE_SOURCE_DIR}/src/algs/Qslim
  ${CMAKE_SOURCE_DIR}/src/algs/Sweep
  ${CMAKE_SOURCE_DIR}/src/extern/CAMAL
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

set( TESTS test_edgemesher.cpp
           test_assygen.cpp
           test_postbl.cpp
           test_meshoptemplate.cpp
           test_facetmeshreader.cpp
           test_facetmeshreader_senses.cpp
           test_tfimapping.cpp
           test_onetooneswept.cpp
           test_ebmesh.cpp
           test_copymesh.cpp
           test_copygeom.cpp
           test_extrudemesh.cpp
           test_scdmesh.cpp
           test_qslimmesher.cpp
    )

## QUADMESHER
# test_quadmeshcleanup.cpp test_jaalquad.cpp
## LPSOLVER
# test_submapping.cpp
## DAGMC
# test_make_watertight.cpp DEPS:src/algs/make_watertight/libMakeWatertight
## FBIGEOM
# test_fbgeom.cpp test_mbgeomop.cpp test_mbsplitop.cpp test_mbvolop.cpp
## CAMAL + Triangle
# test_setpnt2quad.cpp
## ENABLE_intassign
# test_intassign.cpp
## HAVE_CAMAL_PAVER
# test_assymesher.cpp

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

#add_subdirectory( AdvFront )

