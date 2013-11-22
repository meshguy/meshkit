#MOAB supports CMake
add_external_project(moab
  DEPENDS hdf5
  USE_AUTOCONF
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
    --prefix=<INSTALL_DIR>
    --with-hdf5=<INSTALL_DIR>
    --with-zlib=<INSTALL_DIR>
    --with-szip=<INSTALL_DIR>
    --without-damsel
    --without-ccmio
    --enable-shared
)