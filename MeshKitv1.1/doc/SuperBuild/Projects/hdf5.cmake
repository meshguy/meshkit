
get_include_dir(szip szip_include_dir)
get_libraries(szip szip_libraries)

#hdf5 supports CMake
add_external_project(hdf5
  DEPENDS zlib szip
  CMAKE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=TRUE
    -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=TRUE
    -DHDF5_ENABLE_SZIP_SUPPORT:BOOL=TRUE
    -DHDF5_ENABLE_SZIP_ENCODING:BOOL=TRUE
    -DHDF5_BUILD_HL_LIB:BOOL=TRUE
    -DSZIP_LIBRARY:FILEPATH=${szip_libraries}
    -DSZIP_INCLUDE_DIR:FILEPATH=${szip_include_dir}
)

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  ExternalProject_Add_Step(hdf5 cpDebugDylib
    COMMAND ${CMAKE_COMMAND} -E create_symlink
    <INSTALL_DIR>/lib/libhdf5_debug${CMAKE_SHARED_LIBRARY_SUFFIX}
    <INSTALL_DIR>/lib/libhdf5${CMAKE_SHARED_LIBRARY_SUFFIX}
    DEPENDEES install
    )
endif()

ExternalProject_Get_Property(${name} install_dir)
set_libraries(hdf5 ${install_dir}/lib/libhdf5${CMAKE_SHARED_LIBRARY_SUFFIX})
