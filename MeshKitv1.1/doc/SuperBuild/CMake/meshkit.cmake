add_external_project(meshkit
  DEPENDS moab cgm lasso
  USE_AUTOCONF
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
  --prefix=<INSTALL_DIR>
  --with-itaps=<INSTALL_DIR>
  --enable-algs
  --enable-optimize
  --enable-src
  --enable-utils
  --enable-rgg
  --enable-shared
)
#  --with-camal=DIR        Specify location of CAMAL library.
#  --with-tetgen=DIR       Specify location of TetGen library.
#  --with-netgen=DIR       Specify location of Netgen library.
#  --with-triangle=DIR     Specify location of Triangle library.
#  --with-mesquite[=DIR]   Specify where to find Mesquite library
