# This maintains the links for all sources used by this superbuild.
# Simply update this file to change the revision.
# One can use different revision on different platforms.
# e.g.
# if (UNIX)
#   ..
# else (APPLE)
#   ..
# endif()


add_revision(mpich2
  URL "http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/1.4.1p1/mpich2-1.4.1p1.tar.gz"
  URL_MD5 b470666749bcb4a0449a072a18e2c204)

add_revision(freetype
  URL "http://paraview.org/files/dependencies/freetype-2.4.8.tar.gz"
  URL_MD5 "5d82aaa9a4abc0ebbd592783208d9c76")

add_revision(ftgl
  URL "http://sourceforge.net/projects/ftgl/files/FTGL%20Source/2.1.3~rc5/ftgl-2.1.3-rc5.tar.gz"
  URL_MD5 fcf4d0567b7de9875d4e99a9f7423633
  )

add_revision(OCE
  GIT_REPOSITORY "https://github.com/robertmaynard/oce.git"
  GIT_TAG "cgm_support"
  )

add_revision(zlib
  URL "http://www.paraview.org/files/v3.12/zlib-1.2.5.tar.gz"
  URL_MD5 c735eab2d659a96e5a594c9e8541ad63)

add_revision(szip
  URL "http://www.hdfgroup.org/ftp/lib-external/szip/2.1/src/szip-2.1.tar.gz"
  URL_MD5 902f831bcefb69c6b635374424acbead)

add_revision(netcdf
  URL "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.0.tar.gz"
  URL_MD5 40c0e53433fc5dc59296ee257ff4a813)

add_revision(hdf5
  URL "http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.8/src/hdf5-1.8.8.tar.gz"
  URL_MD5 1196e668f5592bfb50d1de162eb16cff)

add_revision(moab
  GIT_REPOSITORY https://bitbucket.org/fathomteam/moab.git
  GIT_TAG "bd52ba12517416f4b6d2162696a41583b73d52ed"
  )

add_revision(cgm
  SVN_REPOSITORY https://svn.mcs.anl.gov/repos/ITAPS/cgm/tags/13.1.1
  SVN_TRUST_CERT 1
  )

add_revision(lasso
  SVN_REPOSITORY https://svn.mcs.anl.gov/repos/ITAPS/Lasso/trunk
  SVN_REVISION "-r6103"
  SVN_TRUST_CERT 1
  )

add_revision(meshkit
  SVN_REPOSITORY https://svn.mcs.anl.gov/repos/fathom/MeshKit/tags/MeshKitv1.1
  SVN_TRUST_CERT 1
  )