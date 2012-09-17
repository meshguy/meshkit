cmake_minimum_required(VERSION 2.8)
project(MeshKitSuperBuild)
set (platform)
if (UNIX)
  set (platform "unix")
elseif (APPLE)
  set (platform "apple")
elseif (WIN32)
  set (platform "win32")
endif()

# Setup CMAKE_MODULE_PATH so that platform specific configurations are processed
# before the generic ones.
set (CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/CMake/${platform}"
  "${CMAKE_CURRENT_SOURCE_DIR}/CMake/common"
  "${CMAKE_CURRENT_SOURCE_DIR}/CMake"
  ${CMAKE_MODULE_PATH}
  )

set (MeshKitSuperBuild_CMAKE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

include (MeshKitModules)
include (${CMAKE_CURRENT_SOURCE_DIR}/versions.cmake)

set (install_location "${CMAKE_CURRENT_BINARY_DIR}/install")
set (download_location "${CMAKE_CURRENT_BINARY_DIR}/downloads")

set (ldflags "$ENV{LDFLAGS} -L${install_location}/lib")
set (cppflags "$ENV{CPPFLAGS} -I${install_location}/include -fPIC")
set (cxxflags "$ENV{CXXFLAGS} -I${install_location}/include -fPIC")
set (cflags "$ENV{CFLAGS} -I${install_location}/include")
set (ld_library_path "${install_location}/lib:$ENV{LD_LIBRARY_PATH}")
set (prefix_path "${install_location}")

string (STRIP "${ldflags}" ldflags)
string (STRIP "${cppflags}" cppflags)
string (STRIP "${cxxflags}" cxxflags)
string (STRIP "${cflags}"  cflags)


#currently remove mpich2 support to make it easier to build for now
#add_project(mpich2)

#we currently don't make anything optional as we are just trying to get
#it working for now
add_project(ftgl REQUIRED)
add_project(OCE REQUIRED)
add_project(zlib REQUIRED)
add_project(szip REQUIRED)
add_project(hdf5 REQUIRED)
add_project(moab REQUIRED)
add_project(cgm REQUIRED)
add_project(lasso REQUIRED)
add_project(meshkit REQUIRED)
