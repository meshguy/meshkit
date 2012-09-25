add_external_project(lasso
  DEPENDS moab cgm
  USE_AUTOCONF
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
    --enable-encoding
    --enable-shared
    --with-imesh=<INSTALL_DIR>
    --with-igeom=<INSTALL_DIR>
    --prefix=<INSTALL_DIR>
)
