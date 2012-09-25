
get_include_dir(OCE OCE_include_dir)

#strip the include directory
get_filename_component(oce_path "${OCE_include_dir}../" PATH)

add_external_project(cgm
  DEPENDS OCE
  USE_AUTOCONF
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
    --with-occ=${oce_path}
    --prefix=<INSTALL_DIR>
    --enable-shared
)
