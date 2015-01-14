set(MakeWatertight_srcs
  arc.cpp
  cleanup.cpp
  gen.cpp
  zip.cpp
  cw_func.cpp
  mw_func.cpp)

set(MakeWatertight_headers
  meshkit/arc.hpp
  meshkit/cleanup.hpp
  meshkit/gen.hpp
  meshkit/zip.hpp
  meshkit/cw_func.hpp
  meshkit/mw_func.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")
if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

include_directories(${MOAB_INCLUDE_DIRS})

add_library(MakeWatertight STATIC
  ${MakeWatertight_srcs}
  ${MakeWatertight_headers})
target_link_libraries(MakeWatertight
  MOAB)

add_definitions(-DTEST)

add_executable(make_watertight
  make_watertight.cpp)
target_link_libraries(make_watertight
  MakeWatertight)

add_executable(check_watertight
  check_watertight.cpp)
target_link_libraries(check_watertight
  MakeWatertight)

install(
  TARGETS     make_watertight
              check_watertight
  DESTINATION bin
  COMPONENT   runtime)
