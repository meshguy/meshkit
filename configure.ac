################################################################################
#                           Standard Stuff
################################################################################
AC_INIT(MESHKIT, 3.99)
AC_CANONICAL_SYSTEM
AC_CONFIG_MACRO_DIR([m4])
m4_ifdef([AM_SILENT_RULES],
 [AM_INIT_AUTOMAKE([-Wall color-tests parallel-tests])],
 [AM_INIT_AUTOMAKE([-Wall])])

dnl Original MK_CHECK_COMPILERS did not look for Fortran.  We
dnl ask FATHOM_CHECK_COMPILERS to also skip the check for FC:
dnl FATHOM_CHECK_COMPILERS([yes],[yes],[no])
dnl However, there is a call to AC_FC_WRAPPERS below, so presumably
dnl we should be checking for FC. - j.kraftcheck 2010-8-13
EXTRA_GNU_FLAGS="-pipe -Wall"
FATHOM_CHECK_COMPILERS
AM_CONDITIONAL(PARALLEL,[test "x$WITH_MPI" != "xno"])
LIBS="-lm"

AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_DISABLE_SHARED
AC_PROG_LIBTOOL
FATHOM_COMPILER_FLAGS
m4_ifdef([AM_SILENT_RULES],[
if test "x$DEBUG" = "xyes"; then
  AM_SILENT_RULES(yes)
else
  AM_SILENT_RULES(no)
fi
])

AC_CHECK_PROG([ZCAT],[gunzip],[gunzip -c],[])
AC_CHECK_PROG([ZCAT],[gzip],[gzip -cd],[])
AC_CHECK_PROG([ZCAT],[zcat],[zcat],[])
if test "x" = "x$ZCAT"; then
  AC_MSG_WARN([[Cannot run tests: no way to uncompress input files.]])
fi
AC_SUBST([ZCAT])
AM_CONDITIONAL([HAVE_ZCAT],[test "x" = "x$ZCAT"])

dnl This doesn't appear to be defined anywhere. - j.kraftcheck 2010-8-13
dnl MK_TEMPLATE_SPECIALIZATION
dnl DEFINES="$DEFINES $TEMPLATE_SPECIALIZATION"

if test "xyes" = "x$ENABLE_FORTRAN" && test "x" != "x$FC"; then
  AC_FC_WRAPPERS
fi

################################################################################
#                           Basic Portability Stuff
################################################################################
AC_CHECK_FUNC([vsnprintf],[AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_VSNPRINTF"])

# Use 64-bit off_t value for file offsets.  W/out this, the stat/fseek calls
# in serial_open_file and other places will fail on 32-bit systems if the file 
# is larger that 2GB.
AM_CPPFLAGS="$AM_CPPFLAGS -D_FILE_OFFSET_BITS=64"



################################################################################
#                              Version Stuff
################################################################################

AC_DEFINE(MK_VERSION,["AC_PACKAGE_VERSION"],[MeshKit Version])
VERSION_MAJOR=`expr AC_PACKAGE_VERSION : '\([[0-9]]*\)'`
VERSION_MINOR=`expr AC_PACKAGE_VERSION : '[[0-9]]*\.\([[0-9]]*\)'`
VERSION_PATCH=`expr AC_PACKAGE_VERSION : '[[0-9]]*\.[[0-9]]*\.\(.*\)'`
test "x" != "x$VERSION_MAJOR" || AC_MSG_ERROR("Invalid version string: AC_PACKAGE_VERSION")
test "x" != "x$VERSION_MINOR" || AC_MSG_ERROR("Invalid version string: AC_PACKAGE_VERSION")
AC_DEFINE_UNQUOTED(MK_VERSION_MAJOR,$VERSION_MAJOR,[MeshKit Major Version])
AC_DEFINE_UNQUOTED(MK_VERSION_MINOR,$VERSION_MINOR,[MeshKit Minor Version])
if test "x" != "x$VERSION_PATCH"; then
  AC_DEFINE_UNQUOTED(MK_VERSION_PATCH,$VERSION_PATCH,[MeshKit Patch Level])
  VERSION_STRING="AC_PACKAGE_STRING"
elif test $VERSION_MINOR -eq 99; then
  VERSION_STRING="AC_PACKAGE_STRING (Alpha)"
else
  VERSION_STRING="AC_PACKAGE_STRING (Beta)"
fi
AC_DEFINE_UNQUOTED(MK_VERSION_STRING,"${VERSION_STRING}",[MeshKit Version String])


################################################################################
#                             System Headers
################################################################################
AC_CHECK_HEADER([ieeefp.h],[DEFINES="$DEFINES -DHAVE_IEEEFP_H"] )
if test "xyes" = "x$enable_debug"; then
  AC_CHECK_HEADER( [valgrind/memcheck.h], [CPPFLAGS="$CPPFLAGS -DVALGRIND"] )
fi


#################################################################################
#                             Doxygen Documentation
#################################################################################

# List of source directories to conditionally search.
# Paths should be relative to top source directory.
DOXYGEN_SRC_PATH="src/core src/utils src/algs src/extern/register_extern.cpp test/core test/algs"

# Old versions of autoconf (automake?) don't define this
if test "x$htmldir" = "x"; then
  htmldir="$docdir"
  AC_SUBST(htmldir)
fi

AC_ARG_ENABLE([html-docs],
[AC_HELP_STRING([[--enable-html-docs]],[Generate and install HTML API documentation])
AC_HELP_STRING([--disable-html-docs],[Do not generate API documentation (default)])],
                        [ENABLE_DOXYGEN="$enableval"],[ENABLE_DOXYGEN=no] )
AC_PATH_PROG( [DOXYGEN], [doxygen], [no] )
if test "x$ENABLE_DOXYGEN" != "xno"; then
  if test "x$DOXYGEN" = "xno"; then
    AC_MSG_ERROR("Doxygen executable not found.")
  fi
fi
AC_ARG_VAR([DOXYGEN],[Doxygen source code documentation tool])
AC_SUBST([DOXYGEN])
AM_CONDITIONAL([ENABLE_DOXYGEN],[test "x$ENABLE_DOXYGEN" != "xno"])


################################################################################
#                              MPI OPTIONS
################################################################################

if test "x$WITH_MPI" != "xno"; then
  DEFINES="$DEFINES -DUSE_MPI"
  AC_CHECK_HEADER([mpi.h],[],[AC_MSG_ERROR([mpi.h not found or not working])])
  AC_MSG_CHECKING([If -DMPICH_IGNORE_CXX_SEEK is necessary])
  AC_LANG_PUSH([C++])
  AC_PREPROC_IFELSE([AC_LANG_PROGRAM( [[#include <stdio.h>
     #include <mpi.h>]],[])],
    [AC_MSG_RESULT([no])],
    [save_CPPFLAGS="$CPPFLAGS"
     CPPFLAGS="$CPPFLAGS -DMPICH_IGNORE_CXX_SEEK"
     AC_PREPROC_IFELSE([AC_LANG_PROGRAM( [[#include <stdio.h>
         #include <mpi.h>]],[])],
         [AC_MSG_RESULT([yes]); AM_CPPFLAGS="$AM_CPPLFLAGS -DMPICH_IGNORE_CXX_SEEK" ],
         [AC_MSG_ERROR([mpi.h is not usable])]),
     CPPFLAGS="$save_CPPFLAGS"])
  AC_LANG_POP([C++])
fi
AM_CONDITIONAL(USE_MPI, [test "xno" != "x$WITH_MPI"])

################################################################################
#                           ITAPS interfaces
################################################################################


SF="$DISTCHECK_CONFIGURE_FLAGS"

AC_ARG_WITH(itaps,
[AC_HELP_STRING([--with-itaps=DIR], 
                [specify same dir for --with-imesh, --with-igeom, --with-irel])],
[ITAPS_DIR="$withval"; SF="$SF --with-itaps=$withval"],[ITAPS_DIR=])

AC_ARG_WITH(igeom,
[AC_HELP_STRING([--with-igeom=DIR],[Directory containing iGeom implementation])],
[IGEOM_DIR=$withval; SF="$SF --with-igeom=$withval"],
[IGEOM_DIR="$ITAPS_DIR"; IGEOM_OKAY_IF_MISSING=yes])

AC_ARG_WITH(imesh,
[AC_HELP_STRING([--with-imesh=DIR],[Directory containing iMesh implementation])],
[IMESH_DIR=$withval; SF="$SF --with-imesh=$withval"],
[IMESH_DIR="$ITAPS_DIR"; IMESH_OKAY_IF_MISSING=yes])

AC_ARG_WITH(irel,
[AC_HELP_STRING([--with-irel=DIR],[Directory containing iRel implementation])],
[IREL_DIR=$withval; SF="$SF --with-irel=$withval"],
[IREL_DIR="$ITAPS_DIR"; IREL_OKAY_IF_MISSING=yes])

DISTCHECK_CONFIGURE_FLAGS="$SF"

#                         Define Macro to check for i*-Defs.inc
# Arguments:
#   file prefix (iGeom, iMesh, etc.)
#   var prefix (IGEOM, IMESH, etc.)
#   path
#   warn if not found (warn if 'yes', error otherwise)
#   warning message
# Does the following
#   appends -DHAVE_Ixxxx = yes/no to AM_CPPFLAGS
#   defines HAVE_Ixxxx automake conditional flag
#   sets Ixxxx_CONFIG_OPTIONS to full path to makefile stub
AC_DEFUN([CHECK_ITAPS_DEFS],[
itaps_pfx="$1"
itaps_dir="$3"
warn="$4"
case "x$itaps_dir" in
  xno|x)
    HAVE_$2=no
    ;;
  xyes)
    AC_MSG_ERROR([Expected agrument for --with-$itaps_pfx])
    ;;
  *)
    if ! test -d "$itaps_dir"; then
      AC_MSG_ERROR([$itaps_dir : not a directory])
    fi
    HAVE_$2=yes
    AC_CHECK_FILE([$itaps_dir/lib/${itaps_pfx}-Defs.inc],
     [$2_CONFIG_OPTIONS="$itaps_dir/lib/${itaps_pfx}-Defs.inc"],
     [AC_CHECK_FILE([$itaps_dir/${itaps_pfx}-Defs.inc],
       [$2_CONFIG_OPTIONS="$itaps_dir/${itaps_pfx}-Defs.inc"],
       [HAVE_$2=no])])
    ;;
esac

if test "xyes" = "x$HAVE_$2"; then
  FATHOM_MAKE_INC_VAR([${$2_CONFIG_OPTIONS}],[$2_INCLUDES],[$2_INCLUDES="$make_val"],
                      [AC_MSG_ERROR([$2_INCLUDES not found in $2_CONFIG_OPTIONS])])
  FATHOM_MAKE_INC_VAR([${$2_CONFIG_OPTIONS}],[$2_LDFLAGS],[$2_LDFLAGS="$make_val"],
                      [AC_MSG_ERROR([$2_LDFLAGS not found in $2_CONFIG_OPTIONS])])
  FATHOM_MAKE_INC_VAR([${$2_CONFIG_OPTIONS}],[$2_LIBS],[$2_LIBS="$make_val"],
                      [AC_MSG_ERROR([$2_LIBS not found in $2_CONFIG_OPTIONS])])
  $2_CONFIG_OPTIONS_INCLUDE="include ${$2_CONFIG_OPTIONS}"
  AC_SUBST($2_CONFIG_OPTIONS_INCLUDE)
  $2_DIR=$itaps_dir
  AC_SUBST($2_DIR)
fi

if test "xyes" = "x$HAVE_$2"; then 
  AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_$2 ${$2_INCLUDES}"
  AM_LDFLAGS="$AM_LDFLAGS ${$2_LDFLAGS}"
  LIBS="${$2_LIBS} $LIBS"
elif test "xyes" = "x$warn"; then
  AC_MSG_WARN([$itaps_pfx not found: $5])
else
  AC_MSG_ERROR([$itaps_pfx not found])
fi
AM_CONDITIONAL( HAVE_$2, [test "xyes" == "x$HAVE_$2"] )
])

CHECK_ITAPS_DEFS( iGeom, IGEOM, $IGEOM_DIR, $IGEOM_OKAY_IF_MISSING, [build may fail] )
CHECK_ITAPS_DEFS( iMesh, IMESH, $IMESH_DIR, $IMESH_OKAY_IF_MISSING, [build may fail] )
CHECK_ITAPS_DEFS( iRel,  IREL,  $IREL_DIR,  $IREL_OKAY_IF_MISSING, [camel disabled] )

################################################################################
#                           MOAB
################################################################################
old_LIBS="$LIBS"
old_LDFLAGS="$LDFLAGS"
old_CPPFLAGS="$CPPFLAGS"
LDFLAGS="$LDFLAGS $IMESH_LDFLAGS"
AC_LANG_PUSH([C++])
AC_CHECK_LIB( [MOAB], [main], [HAVE_MOAB=yes],[HAVE_MOAB=no] )
AM_CONDITIONAL( HAVE_MOAB, [test "xyes" == "x$HAVE_$2"] )
AC_SUBST([HAVE_MOAB])
if test "xyes" = "x$HAVE_MOAB"; then 
  AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_MOAB"
fi
LIBS="$old_LIBS"
LDFLAGS="$old_LDFLAGS"
CPPFLAGS="$old_CPPFLAGS"

################################################################################
#                           CGM
################################################################################
old_LIBS="$LIBS"
old_LDFLAGS="$LDFLAGS"
old_CPPFLAGS="$CPPFLAGS"
LDFLAGS="$LDFLAGS $IGEOM_LDFLAGS"
AC_LANG_PUSH([C++])
AC_CHECK_LIB( [cgm], [main], [HAVE_CGM=yes],[HAVE_CGM=no] )
AM_CONDITIONAL( HAVE_CGM, [test "xyes" == "x$HAVE_$2"] )
AC_SUBST([HAVE_CGM])
if test "xyes" = "x$HAVE_CGM"; then 
  AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_CGM"
fi
LIBS="$old_LIBS"
LDFLAGS="$old_LDFLAGS"
CPPFLAGS="$old_CPPFLAGS"

###############################################################################
#                           Optional Tools
###############################################################################
  # Define a macro to avoid typing this for each individual tool
  # Usage: MK_OPTIONAL_TOOL( name, default )
  #  name      - name of option
  #  default   - yes/no, $ENABLE_TOOLS overrides this if set
  # Actions:
  #  sets ENABLE_${tool} to 'yes' or 'no'
  #  creates ENABLE_${tool} automake conditional
AC_DEFUN([MK_OPTIONAL_TOOL],[
  mk_default_$1=$2
  if test "x" != "x$ENABLE_TOOLS"; then
    mk_default_$1="$ENABLE_TOOLS"
  fi
  AC_ARG_ENABLE( [$1],
[AC_HELP_STRING([--enable-$1],[Build tool: $1])
AC_HELP_STRING([--disable-$1],[Don't build $1])],
                 [ENABLE_$1=$enableval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-$1=\"${enableval}\""
],[ENABLE_$1=${mk_default_$1}] )
  AM_CONDITIONAL([ENABLE_$1],[test "x${ENABLE_$1}" != "xno"])
  AC_MSG_CHECKING([if $1 is to be built])
  AC_MSG_RESULT([${ENABLE_$1}])
])
  # The default for all tools
AC_ARG_ENABLE( [tools], [
AC_HELP_STRING([--enable-tools],[Build all tools by default])
AC_HELP_STRING([--disable-tools],[Disable all tools by default])],
                         [ENABLE_TOOLS=$enableval],[ENABLE_TOOLS=] )
  # Individual tools
MK_OPTIONAL_TOOL([algs],   [yes])
MK_OPTIONAL_TOOL([rgg],    [yes])
MK_OPTIONAL_TOOL([utils],  [yes])
MK_OPTIONAL_TOOL([src],    [no])
MK_OPTIONAL_TOOL([python], [no])


################################################################################
#                           CAMAL
################################################################################
AC_ARG_WITH(camal, 
[AC_HELP_STRING([--with-camal=DIR], [Specify location of CAMAL library.])],
[CAMAL_VAL=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-camal=\"${withval}\""
 ], [CAMAL_VAL=no])
 
CAMAL_VERSION=
CAMAL_CPPFLAGS=
CAMAL_LDFLAGS=
CAMAL_LIBS=

# Hard-code versions because CAMAL lib names are broken
# (lib names should not include a version number unless 
# a change in that number will break all code using the
# library, which certainly shouldn't be the case for the
# minor version number or the patch level, and for the
# major version, the standard library versioning scheme
# should be used.)  Further, there is no way to determine 
# the version from the header files, only by linking with
# the library and calling API functions, which requires
# knowing the version.
# NOTE: This is a space-separated list.  Other supported 
#       versions may be appended.
CAMAL_VERSIONS="510 330 310"

if test "xyes" != "x$HAVE_IREL"; then
  case "x$CAMAL_VAL" in
    xno)
      ;;
    x|xyes)
      AC_MSG_WARN([Camal support disabled because iRel not found])
      CAMAL_VAL=no
      ;;
    *)
      AC_MSG_ERROR([Camal cannot be enabled w/out iRel])
      ;;
  esac
fi

if test "xno" = "x$CAMAL_VAL"; then
  HAVE_CAMAL=no
  HAVE_CAMAL_EDGE=no
  HAVE_CAMAL_PAVER=no
  HAVE_CAMAL_SURFMAPPER=no
  HAVE_CAMAL_TRIADV=no
  HAVE_CAMAL_TRIDEL=no
  HAVE_CAMAL_SWEEP=no
  HAVE_CAMAL_TET=no
else
    
  CAMAL_DIR="$CAMAL_VAL"
  # if user didn't specify anything, assume 'yes' for now
  # such that we enable it if we find camal in the default
  # search path (or user specifed CPPFLAGS and LDFLAGS)
  if test "x" = "x$CAMAL_VAL"; then
    CAMAL_DIR=yes
  fi
  
  # if user specified a path, add subdirs to search path
  if test "xyes" != "x$CAMAL_DIR"; then
    if ! test -d "$CAMAL_DIR"; then
      AC_MSG_ERROR([$CAMAL_DIR: not a directory])
    fi
    
    if test -d "$CAMAL_DIR/include"; then
      CAMAL_CPPFLAGS="-I$CAMAL_DIR/include"
    else
      AC_MSG_WARN([$CAMAL_DIR/include: no such directory])
      CAMAL_CPPFLAGS="-I$CAMAL_DIR"
    fi
    
    if test -d "$CAMAL_DIR/lib"; then
      # CAMAL puts libs in subdirs by a rather imprecise mishmash 
      # of OS and arch.  Just add them all to LDFLAGS.  The linker
      # will ignore libraries for other platforms.
      for file in "$CAMAL_DIR"/lib/*; do
          # Skip anything that isn't a directory or that
          # begins with a '.'.
        case $file in
          .*)
            ;;
          *)
            if test -d $file; then
              CAMAL_LDFLAGS="$CAMAL_LDFLAGS -L$file"
            fi
            ;;
        esac
      done
    else
      AC_MSG_WARN([$CAMAL_DIR/lib: no such directory])
      CAMAL_LDFLAGS="-L$CAMAL_DIR"
    fi
  fi
  
  # Check for CAMAL headers and libraries
  old_LIBS="$LIBS"
  old_LDFLAGS="$LDFLAGS"
  old_CPPFLAGS="$CPPFLAGS"
  LDFLAGS="$LDFLAGS $CAMAL_LDFLAGS"
  CPPFLAGS="$CPPFLAGS $CAMAL_CPPFLAGS"
  AC_LANG_PUSH([C++])
  
  
  
  for ver in $CAMAL_VERSIONS; do
    AC_CHECK_LIB( [CMLSizer], [main], [CAMAL_VERSION=$ver; break;] )
    AC_CHECK_LIB( [CMLShr$ver], [main], [CAMAL_VERSION=$ver; break;] )
  done
  if test "x$CAMAL_VERSION" = "x"; then
    HAVE_CAMAL=no
    HAVE_CAMAL_TRIADV=no
    HAVE_CAMAL_PAVER=no
    HAVE_CAMAL_SURFMAPPER = no
    HAVE_CAMAL_SWEEP=no
    HAVE_CAMAL_TET=no
    HAVE_CAMAL_EDGE=no
  else
    if test "x$CAMAL_VERSION" = "x510"; then
      AC_CHECK_HEADER( [CMLTriAdvance.hpp],[HAVE_CAMAL_TRIADV=yes], [HAVE_CAMAL_TRIADV=no] )
      AC_CHECK_HEADER( [CMLTriDelaunay.hpp],[HAVE_CAMAL_TRIDEL=yes], [HAVE_CAMAL_TRIDEL=no] )
      AC_CHECK_HEADER( [CMLPaver.hpp],    [HAVE_CAMAL_PAVER=yes],  [HAVE_CAMAL_PAVER=no] )
      AC_CHECK_HEADER( [CMLSurfMapper.hpp],    [HAVE_CAMAL_SURFMAPPER=yes],  [HAVE_CAMAL_SURFMAPPER=no] )
      AC_CHECK_HEADER( [CMLSweeper.hpp],  [HAVE_CAMAL_SWEEP=yes],  [HAVE_CAMAL_SWEEP=no] )
      AC_CHECK_HEADER( [CMLTetMesher.hpp],[HAVE_CAMAL_TET=yes],    [HAVE_CAMAL_TET=no] )
      AC_CHECK_HEADER( [CMLSurfEval.hpp], [HAVE_CAMAL=yes],        [HAVE_CAMAL=no] )
      AC_CHECK_HEADER( [CMLEdgeMesher.hpp], [HAVE_CAMAL=yes],        [HAVE_CAMAL=no] )

      CAMAL_LIBS="-lCMLSizer"
      AC_CHECK_LIB( [CMLPave],   [main], [CAMAL_LIBS="-lCMLPave $CAMAL_LIBS"], [HAVE_CAMAL_PAVER=no],  [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLSurfMapper],   [main], [CAMAL_LIBS="-lCMLSurfMapper $CAMAL_LIBS"], [HAVE_CAMAL_SURFMAPPER=no],  [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLTri],    [main], [CAMAL_LIBS="-lCMLTri $CAMAL_LIBS"], [HAVE_CAMAL_TRIADV=no], [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLDelTri], [main], [CAMAL_LIBS="-lCMLDelTri $CAMAL_LIBS"], [HAVE_CAMAL_TRIDEL=no], [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLSweep],  [main], [CAMAL_LIBS="-lCMLSweep $CAMAL_LIBS"], [HAVE_CAMAL_SWEEP=no],  [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLTet],    [main], [CAMAL_LIBS="-lCMLTet $CAMAL_LIBS"], [HAVE_CAMAL_TET=no],    [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLEdgeMesher],    [main], [CAMAL_LIBS="-lCMLEdgeMesher $CAMAL_LIBS"], [HAVE_CAMAL_EDGE=no],    [$CAMAL_LIBS] )
    else
      AC_CHECK_HEADER( [CMLTriMesher.hpp],[HAVE_CAMAL_TRIADV=yes], [HAVE_CAMAL_TRIADV=no] )
      AC_CHECK_HEADER( [CMLPaver.hpp],    [HAVE_CAMAL_PAVER=yes],  [HAVE_CAMAL_PAVER=no] )
      AC_CHECK_HEADER( [CMLSurfMapper.hpp],    [HAVE_CAMAL_SURFMAPPER=yes],  [HAVE_CAMAL_SURFMAPPER=no] )
      AC_CHECK_HEADER( [CMLSweeper.hpp],  [HAVE_CAMAL_SWEEP=yes],  [HAVE_CAMAL_SWEEP=no] )
      AC_CHECK_HEADER( [CMLTetMesher.hpp],[HAVE_CAMAL_TET=yes],    [HAVE_CAMAL_TET=no] )
      AC_CHECK_HEADER( [CMLGeomEval.hpp], [HAVE_CAMAL=yes],        [HAVE_CAMAL=no] )

      CAMAL_LIBS="-lCMLShr$CAMAL_VERSION -lCMLUtl$CAMAL_VERSION -lCMLShr$CAMAL_VERSION"
      AC_CHECK_LIB( [CMLPave$CAMAL_VERSION],  [main], [CAMAL_LIBS="-lCMLPave$CAMAL_VERSION  $CAMAL_LIBS"], [HAVE_CAMAL_PAVER=no],  [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLSurfMapper$CAMAL_VERSION],   [main], [CAMAL_LIBS="-lCMLSurfMapper$CAMAL_VERSION $CAMAL_LIBS"], [HAVE_CAMAL_SURFMAPPER=no],  [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLTri$CAMAL_VERSION],   [main], [CAMAL_LIBS="-lCMLTri$CAMAL_VERSION   $CAMAL_LIBS"], [HAVE_CAMAL_TRIADV=no], [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLSweep$CAMAL_VERSION], [main], [CAMAL_LIBS="-lCMLSweep$CAMAL_VERSION $CAMAL_LIBS"], [HAVE_CAMAL_SWEEP=no],  [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLTet$CAMAL_VERSION],   [main], [CAMAL_LIBS="-lCMLTet$CAMAL_VERSION   $CAMAL_LIBS"], [HAVE_CAMAL_TET=no],    [$CAMAL_LIBS] )
    fi
    
      # CMLTet requires libg2c of built with g77: check if it exists on this system
    AC_CHECK_LIB( [g2c], [s_cmp], [CAMAL_LIBS="$CAMAL_LIBS -lg2c"] )
    
  fi
  
  AC_LANG_POP([C++])
  LIBS="$old_LIBS"
  LDFLAGS="$old_LDFLAGS"
  CPPFLAGS="$old_CPPFLAGS"
fi

if test $HAVE_CAMAL = yes; then
  AC_MSG_RESULT([Camal version: $CAMAL_VERSION])
  DOXYGEN_SRC_PATH="$DOXYGEN_SRC_PATH src/extern/CAMAL test/extern/CAMAL"
elif test "x$CAMAL_VAL" != "xno"; then
  AC_MSG_ERROR([Camal not found or not usable.])
fi

AC_SUBST([CAMAL_VERSION])
AC_SUBST([CAMAL_CPPFLAGS])
AC_SUBST([CAMAL_LDFLAGS])
AC_SUBST([CAMAL_LIBS])

AC_SUBST([HAVE_CAMAL])
AC_SUBST([HAVE_CAMAL_TRIADV])
AC_SUBST([HAVE_CAMAL_TRIDEL])
AC_SUBST([HAVE_CAMAL_PAVER])
AC_SUBST([HAVE_CAMAL_SWEEP])
AC_SUBST([HAVE_CAMAL_TET])
AC_SUBST([HAVE_CAMAL_EDGE])

AM_CONDITIONAL( HAVE_CAMAL,        [test "x$HAVE_CAMAL"        = "xyes"] )
AM_CONDITIONAL( HAVE_CAMAL_TRIADV, [test "x$HAVE_CAMAL_TRIADV" = "xyes"] )
AM_CONDITIONAL( HAVE_CAMAL_TRIDEL, [test "x$HAVE_CAMAL_TRIDEL" = "xyes"] )
AM_CONDITIONAL( HAVE_CAMAL_PAVER,  [test "x$HAVE_CAMAL_PAVER"  = "xyes"] )
AM_CONDITIONAL( HAVE_CAMAL_SWEEP,  [test "x$HAVE_CAMAL_SWEEP"  = "xyes"] )
AM_CONDITIONAL( HAVE_CAMAL_TET,    [test "x$HAVE_CAMAL_TET"    = "xyes"] )
AM_CONDITIONAL( HAVE_CAMAL_EDGE,   [test "x$HAVE_CAMAL_EDGE"   = "xyes"] )

if test "xyes" = "x$HAVE_CAMAL"; then 
  AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_CAMAL"
  if test "xyes" = "x$HAVE_CAMAL_TRIADV"; then 
    AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_CAMAL_TRIADV"
  fi
  if test "xyes" = "x$HAVE_CAMAL_TRIDEL"; then 
    AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_CAMAL_TRIDEL"
  fi
  if test "xyes" = "x$HAVE_CAMAL_PAVER"; then 
    AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_CAMAL_PAVER"
  fi
  if test "xyes" = "x$HAVE_CAMAL_SWEEP"; then 
    AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_CAMAL_SWEEP"
  fi
  if test "xyes" = "x$HAVE_CAMAL_TET"; then 
    AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_CAMAL_TET"
  fi
  if test "xyes" = "x$HAVE_CAMAL_EDGE"; then 
    AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_CAMAL_EDGE"
  fi
fi

################################################################################
#                           TetGen
################################################################################
AC_ARG_WITH(tetgen, 
[AC_HELP_STRING([--with-tetgen=DIR], [Specify location of TetGen library.])],
[TETGEN_VAL=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-tetgen=\"${withval}\""
 ], [TETGEN_VAL=no])
 
TETGEN_CPPFLAGS=
TETGEN_LDFLAGS=
TETGEN_LIBS=

# check if iRel is enable
if test "xyes" != "x$HAVE_IREL"; then
  case "x$TETGEN_VAL" in
    xno)
      ;;
    x|xyes)
      AC_MSG_WARN([TetGen support disabled because iRel not found])
      TETGEN_VAL=no
      ;;
    *)
      AC_MSG_ERROR([TetGen cannot be enabled w/out iRel])
      ;;
  esac
fi

# Check for TETGEN headers and libraries
if test "xno" = "x$TETGEN_VAL"; then
  HAVE_TETGEN=no
else
  TETGEN_DIR="$TETGEN_VAL"

  if test "x" = "x$TETGEN_VAL"; then
    TETGEN_DIR=yes
  fi

  if test "xyes" != "x$TETGEN_DIR"; then
    if ! test -d "$TETGEN_DIR"; then
      AC_MSG_ERROR([$TETGEN_DIR: not a directory])
    fi
    
    if test -d "$TETGEN_DIR/include"; then
      TETGEN_CPPFLAGS="-I$TETGEN_DIR/include"
    else
      AC_MSG_WARN([$TETGEN_DIR/include: no such directory])
      TETGEN_CPPFLAGS="-I$TETGEN_DIR"
    fi
    
    if test -d "$TETGEN_DIR/lib"; then
      TETGEN_LDFLAGS="-L$TETGEN_DIR/lib"
    else
      AC_MSG_WARN([$TETGEN_DIR/lib: no such directory])
      TETGEN_LDFLAGS="-L$TETGEN_DIR"
    fi
  fi

  old_LIBS="$LIBS"
  old_LDFLAGS="$LDFLAGS"
  old_CPPFLAGS="$CPPFLAGS"
  LDFLAGS="$LDFLAGS $TETGEN_LDFLAGS -ltet"
  CPPFLAGS="$CPPFLAGS $TETGEN_CPPFLAGS"
  AC_LANG_PUSH([C++])

  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include "tetgen.h"],
        [tetgenio dummy])],
    [TETGEN_LIBS="$TETGEN_LDFLAGS -ltet"] [HAVE_TETGEN=yes],
    [AC_MSG_WARN([libtet is not installed.])])

  AC_LANG_POP([C++])
  LIBS="$old_LIBS"
  LDFLAGS="$old_LDFLAGS"
  CPPFLAGS="$old_CPPFLAGS"
fi

if test "x$HAVE_TETGEN" = "xyes"; then
  AC_MSG_RESULT([TetGen is found.])
elif test "x$TETGEN_VAL" != "xno"; then
  AC_MSG_ERROR([TetGen not found or not usable.])
fi

AC_SUBST([TETGEN_CPPFLAGS])
AC_SUBST([TETGEN_LDFLAGS])
AC_SUBST([TETGEN_LIBS])
AC_SUBST([HAVE_TETGEN])

AM_CONDITIONAL(HAVE_TETGEN, [test "x$HAVE_TETGEN" = "xyes"])


################################################################################
#                           NetGen
################################################################################
AC_ARG_WITH(netgen, 
[AC_HELP_STRING([--with-netgen=DIR], [Specify location of Netgen library.])],
[NETGEN_VAL=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-netgen=\"${withval}\""
 ], [NETGEN_VAL=no])
 
NETGEN_CPPFLAGS=
NETGEN_LDFLAGS=
NETGEN_LIBS=

# Check for NETGEN headers and libraries
if test "xno" = "x$NETGEN_VAL"; then
  HAVE_NETGEN=no
else
  NETGEN_DIR="$NETGEN_VAL"

  if test "x" = "x$NETGEN_VAL"; then
    NETGEN_DIR=yes
  fi

  if test "xyes" != "x$NETGEN_DIR"; then
    if ! test -d "$NETGEN_DIR"; then
      AC_MSG_ERROR([$NETGEN_DIR: not a directory])
    fi
    
    if test -d "$NETGEN_DIR/include"; then
      NETGEN_CPPFLAGS="-I$NETGEN_DIR/include"
    else
      AC_MSG_WARN([$NETGEN_DIR/include: no such directory])
      NETGEN_CPPFLAGS="-I$NETGEN_DIR"
    fi
    
    if test -d "$NETGEN_DIR/lib"; then
      NETGEN_LDFLAGS="-L$NETGEN_DIR/lib"
    else
      AC_MSG_WARN([$NETGEN_DIR/lib: no such directory])
      NETGEN_LDFLAGS="-L$NETGEN_DIR"
    fi
  fi

  old_LIBS="$LIBS"
  old_LDFLAGS="$LDFLAGS"
  old_CPPFLAGS="$CPPFLAGS"
  LDFLAGS="$LDFLAGS $NETGEN_LDFLAGS -lnglib"
  CPPFLAGS="$CPPFLAGS $NETGEN_CPPFLAGS"
  AC_LANG_PUSH([C++])

  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include "nglib.h"],
        [Ng_Mesh mesh])],
    [NETGEN_LIBS="$NETGEN_LDFLAGS -lnglib"] [HAVE_NETGEN=yes],
    [AC_MSG_WARN([libnglib is not installed.])])

  AC_LANG_POP([C++])
  LIBS="$old_LIBS"
  LDFLAGS="$old_LDFLAGS"
  CPPFLAGS="$old_CPPFLAGS"
fi

if test "x$HAVE_NETGEN" = "xyes"; then
  AC_MSG_RESULT([Netgen is found.])
elif test "x$NETGEN_VAL" != "xno"; then
  AC_MSG_ERROR([Netgen not found or not usable.])
fi

AC_SUBST([NETGEN_CPPFLAGS])
AC_SUBST([NETGEN_LDFLAGS])
AC_SUBST([NETGEN_LIBS])
AC_SUBST([HAVE_NETGEN])

AM_CONDITIONAL(HAVE_NETGEN, [test "x$HAVE_NETGEN" = "xyes"])


################################################################################
#                           Triangle
################################################################################
AC_ARG_WITH(triangle, 
[AC_HELP_STRING([--with-triangle=DIR], [Specify location of Triangle library.])],
[TRIANGLE_VAL=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-triangle=\"${withval}\""
 ], [TRIANGLE_VAL=no])

TRIANGLE_CPPFLAGS=
TRIANGLE_LIBS=

# Check for TRIANGLE headers and libraries
if test "xno" = "x$TRIANGLE_VAL"; then
  HAVE_TRIANGLE=no
else
  TRIANGLE_DIR="$TRIANGLE_VAL"

  if test "x" = "x$TRIANGLE_VAL"; then
    TRIANGLE_DIR=yes
  fi

  if test "xyes" != "x$TRIANGLE_DIR"; then
    if ! test -d "$TRIANGLE_DIR"; then
      AC_MSG_ERROR([$TRIANGLE_DIR: not a directory])
    fi
    TRIANGLE_CPPFLAGS="-I$TRIANGLE_DIR" 
    
  fi
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $TRIANGLE_CPPFLAGS"
  AC_LANG_PUSH([C++])

  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([
        [#define REAL double]
        [#include "triangle.h"]
        ],
        [triangulateio dummy])],
    [TRIANGLE_LIBS="$TRIANGLE_DIR/triangle.o"] [HAVE_TRIANGLE=yes],
    [AC_MSG_WARN([triangle.o is not installed.])])

  AC_LANG_POP([C++])
  LDFLAGS="$old_LDFLAGS"
  CPPFLAGS="$old_CPPFLAGS"

fi

if test "x$HAVE_TRIANGLE" = "xyes"; then
  AC_MSG_RESULT([Triangle is found.])
elif test "x$TRIANGLE_VAL" != "xno"; then
  AC_MSG_ERROR([Triangle not found or not usable.])
fi

AC_SUBST([TRIANGLE_CPPFLAGS])
AC_SUBST([TRIANGLE_LIBS])
AC_SUBST([HAVE_TRIANGLE])

AM_CONDITIONAL(HAVE_TRIANGLE, [test "x$HAVE_TRIANGLE" = "xyes"])


################################################################################
#                           Mesquite
################################################################################

AC_ARG_WITH([mesquite],
            [AC_HELP_STRING([--with-mesquite@<:@=DIR@:>@], [Specify where to find Mesquite library])
AC_HELP_STRING([--without-mesquite], [Configure w/out support for Mesquite library])],
            [WITH_MESQUITE=$withval],
            [WITH_MESQUITE=])


# define function that tries to find a no-arg constructor for a 
# mesquite object in specified lib.  Arguments are 1) class name, 2) library name,
# 3) actions if success and 4) actions if not found.
AC_DEFUN([MESQUITE_LINK_DESTRUCTOR],[
    msq_link_dest_old_LIBS="$LIBS"
    LIBS="-l$2 $LIBS"
    AC_MSG_CHECKING([for $1 in -l$2])
    AC_LINK_IFELSE( [AC_LANG_PROGRAM(
                      [#include "Mesquite.hpp"
                       namespace MESQUITE_NS {
                       class $1 { public: ~$1 (); };
                       }],
                      [MESQUITE_NS::$1 *i; i->~$1()])],
                     [AC_MSG_RESULT([yes])
                      LIBS="$msq_link_dest_old_LIBS"
                      $3],
                     [AC_MSG_RESULT([no])
                      LIBS="$msq_link_dest_old_LIBS"
                      $4] )
])


if test "x$WITH_MESQUITE" != "xno"; then
  AC_LANG_PUSH([C++])

  MESQUITE_LDFLAGS=
  MESQUITE_CPPFLAGS=
  MESQUITE_LIBS=
  HAVE_MESQUITE=yes
  MESQUITE_HAS_IMESH=no
  MESQUITE_HAS_IGEOM=no
  MESQUITE_HAS_IREL=no
  if test "x$WITH_MESQUITE" != "xyes"; then
    if test "x$WITH_MESQUITE" != "x"; then
      MESQUITE_LDFLAGS=-L$WITH_MESQUITE/lib
      MESQUITE_CPPFLAGS=-I$WITH_MESQUITE/include
    fi
  fi
  
  old_LIBS="$LIBS"
  old_LDFLAGS="$LD_FLAGS"
  old_CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$LDFLAGS $MESQUITE_LDFLAGS"
  CPPFLAGS="$CPPFLAGS $MESQUITE_CPPFLAGS"
  
  AC_MSG_CHECKING([for Mesquite.hpp])
  AC_CHECK_HEADER([Mesquite.hpp], ,[HAVE_MESQUITE=no])
  
  if test $HAVE_MESQUITE = yes; then
    MESQUITE_LINK_DESTRUCTOR([PatchData],[mesquite],[LIBS="-lmesquite $LIBS"
                            MESQUITE_LIBS="-lmesquite"],[HAVE_MESQUITE=no])
  fi
  if test $HAVE_MESQUITE = yes; then
    
    MESQUITE_LINK_DESTRUCTOR([MsqIMesh],[mesquite],[MESQUITE_HAS_IMESH=yes],
      [MESQUITE_LINK_DESTRUCTOR([MsqIMesh],[msqimesh], 
        [MESQUITE_HAS_IMESH=yes 
         MESQUITE_LIBS="-lmsqimesh $MESQUITE_LIBS"
         LIBS="-lmsqimesh $LIBS"],[MESQUITE_HAS_IMESH=no])])
    MESQUITE_LINK_DESTRUCTOR([MsqIGeom],[mesquite],[MESQUITE_HAS_IGEOM=yes],
      [MESQUITE_LINK_DESTRUCTOR([MsqIGeom],[msqigeom], 
        [MESQUITE_HAS_IGEOM=yes 
         MESQUITE_LIBS="-lmsqigeom $MESQUITE_LIBS"
         LIBS="-lmsqigeom $LIBS"],[MESQUITE_HAS_IGEOM=no])])
    MESQUITE_LINK_DESTRUCTOR([MsqIRel],[mesquite],[MESQUITE_HAS_IREL=yes],
      [MESQUITE_LINK_DESTRUCTOR([MsqIRel],[msqirel], 
        [MESQUITE_HAS_IREL=yes 
         MESQUITE_LIBS="-lmsqirel $MESQUITE_LIBS"
         LIBS="-lmsqirel $LIBS"],[MESQUITE_HAS_IREL=no])])
  
  fi
  
  LIBS="$old_LIBS"
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  AC_LANG_POP([C++])
  
  if test "xyes" = "x$HAVE_MESQUITE"; then
    if test "xno" = "x$MESQUITE_HAS_IMESH"; then
      AC_MSG_WARN([Mesquite library does not have iMesh support])
      HAVE_MESQUITE=no
    fi
  fi
  
  if test "xno" = "x$HAVE_MESQUITE"; then
    if test "x" != "x$WITH_MESQUITE"; then
      AC_MSG_ERROR([Mesquite not found or not usable])
    fi
    DOXYGEN_SRC_PATH="$DOXYGEN_SRC_PATH src/extern/mesquite test/extern/mesquite"
  fi

  MESQUITE_DEFINES=
  test $MESQUITE_HAS_IMESH = yes && MESQUITE_DEFINES="$MESQUITE_DEFINES -DMSQIMESH"
  test $MESQUITE_HAS_IGEOM = yes && MESQUITE_DEFINES="$MESQUITE_DEFINES -DMSQIGEOM"
  test $MESQUITE_HAS_IREL  = yes && MESQUITE_DEFINES="$MESQUITE_DEFINES -DMSQIREL"
  AC_SUBST([MESQUITE_DEFINES])
  AC_SUBST([MESQUITE_LDFLAGS])
  AC_SUBST([MESQUITE_CPPFLAGS])
  AC_SUBST([MESQUITE_LIBS])
  AM_CONDITIONAL(HAVE_MESQUITE, [test $HAVE_MESQUITE = yes])
fi



################################################################################
#                           Output Files
################################################################################

# Build complete doxygen source path list
TMP_DOXY_LIST=`echo $DOXYGEN_SRC_PATH | $SED -e "s|[[^[:space:]]]*|../&|g"`
TMP_DOXY_LIST2=
if test "x$srcdir" != "x."; then
  TMP_DOXY_LIST2=`echo $DOXYGEN_SRC_PATH | $SED -e "s|[[^[:space:]]]*|../${srcdir}/&|g"`
fi
DOXYGEN_SRC_PATH="$TMP_DOXY_LIST $TMP_DOXY_LIST2"
AC_SUBST(DOXYGEN_SRC_PATH)

AC_SUBST([DEFINES])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([DISTCHECK_CONFIGURE_FLAGS])

AC_ARG_VAR([FC], [FORTRAN compiler command])
AC_CONFIG_HEADERS([MKVersion.h])
AC_CONFIG_FILES([Makefile 
                 meshkit.make 
                 testdir.h
                 algs/Makefile
                 algs/test/Makefile
                 algs/Qslim/Makefile 
                 algs/make_watertight/Makefile
                 camal/Makefile
                 plibs/Makefile
                 python/setup.py
                 rgg/Makefile
                 src/Makefile
                 src/core/Makefile
                 src/algs/Makefile
                 src/algs/QuadMesher/Makefile
                 src/algs/Qslim/Makefile
                 src/utils/Makefile
                 src/lemon/Makefile
                 src/lemon/doc/Makefile
                 src/lemon/doc/Doxyfile
                 src/extern/Makefile
                 src/extern/CAMAL/Makefile
                 src/extern/netgen/Makefile
                 src/extern/mesquite/Makefile
                 test/Makefile
                 test/core/Makefile
                 test/algs/Makefile
                 test/extern/Makefile
                 test/extern/CAMAL/Makefile
                 test/extern/mesquite/Makefile
                 test/extern/netgen/Makefile
                 utils/Makefile
                 doc/Makefile
                 doc/meshkit.dox
                 ])

AC_OUTPUT

AC_MSG_RESULT([C:   $CC $CFLAGS])
AC_MSG_RESULT([C++: $CXX $CXXFLAGS])
AC_MSG_RESULT([Fortran: $FC $FCFLAGS])

if test "x$IMESH_MISSING" = "xyes"; then
  AC_MSG_WARN([Could not find iMesh library.  Build will FAIL!!! Try \"--with-imesh=DIR\"])
fi
