set(Qslim_srcs
  3D.cpp
  Heap.hpp
  Mat4.h
  QslimDecimation.hpp
  std.h
  3D.h
  QslimMesher.cpp
  Vec3.h
  Array.h
  QslimOptions.cpp
  Vec4.h
  Buffer.h
  primitives.cpp
  defs.h
  primitives.h
  quadrics.cpp
  Heap.cpp
  Mat4.cpp
  QslimDecimation.cpp
  quadrics.h)

set(Qslim_headers
  meshkit/QslimMesher.hpp
  meshkit/QslimOptions.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")

add_library(Qslim
  ${Qslim_srcs}
  ${Qslim_headers})

link_libraries(Qslim
  MeshKitcore)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/lemon
    )

install(
  TARGETS   "Qslim"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
