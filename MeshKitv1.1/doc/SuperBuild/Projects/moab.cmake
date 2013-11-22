#MOAB supports CMake

get_libraries(szip szip_libraries)
get_libraries(zlib zlib_libraries)
get_libraries(hdf5 hdf5_libraries)

#strip the file name and the lib directory for moab
get_filename_component(szip_path ${szip_libraries} PATH)
get_filename_component(zlib_path ${zlib_libraries} PATH)
get_filename_component(hdf5_path ${hdf5_libraries} PATH)

#strip the lib directory
get_filename_component(szip_path "${szip_path}../" PATH)
get_filename_component(zlib_path "${zlib_path}../" PATH)
get_filename_component(hdf5_path "${hdf5_path}../" PATH)

#we are presuming the library for all these projects are
#in a folder called lib, on somethind systems like
#ubuntu that isn't true, so moab breaks.
add_external_project(moab
  DEPENDS hdf5
  USE_AUTOCONF
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
    --prefix=<INSTALL_DIR>
    --with-hdf5=${hdf5_path}
    --with-zlib=${zlib_path}
    --with-szip=${szip_path}
    --without-damsel
    --without-ccmio
    --enable-shared
)