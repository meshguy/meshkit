
add_external_project(cgm
  DEPENDS OCE
  USE_AUTOCONF
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
    --with-occ=<INSTALL_DIR>
    --prefix=<INSTALL_DIR>
    --enable-shared
)
