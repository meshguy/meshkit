add_external_project(szip
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
    --enable-encoding
    --prefix=<INSTALL_DIR>
)
