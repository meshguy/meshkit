add_external_project(szip
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
    --enable-encoding
    --prefix=<INSTALL_DIR>
)

ExternalProject_Get_Property(${name} install_dir)
set_libraries(szip ${install_dir}/lib/libsz${CMAKE_SHARED_LIBRARY_SUFFIX})
