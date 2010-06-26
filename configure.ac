################################################################################
#                           Standard Stuff
################################################################################
AC_INIT(MESHKIT, 3.99)
AC_CANONICAL_SYSTEM
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

MK_CHECK_COMPILERS
AM_CONDITIONAL(PARALLEL,[test "x$WITH_MPI" != "xno"])
LIBS="-lm"

AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

AC_CHECK_PROG([ZCAT],[gunzip],[gunzip -c],[])
AC_CHECK_PROG([ZCAT],[gzip],[gzip -cd],[])
AC_CHECK_PROG([ZCAT],[zcat],[zcat],[])
if test "x" = "x$ZCAT"; then
  AC_MSG_WARN([[Cannot run tests: no way to uncompress input files.]])
fi
AC_SUBST([ZCAT])
AM_CONDITIONAL([HAVE_ZCAT],[test "x" = "x$ZCAT"])

MK_TEMPLATE_SPECIALIZATION
DEFINES="$DEFINES $TEMPLATE_SPECIALIZATION"

if test "x" != "x$FC"; then
  AC_FC_WRAPPERS
fi

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

################################################################################
#                              MPI OPTIONS
################################################################################

if test "x$WITH_MPI" != "xno"; then
  DEFINES="$DEFINES -DUSE_MPI"
  AC_CHECK_HEADER([mpi.h],[],[AC_MSG_ERROR([mpi.h not found or not working])])
  AC_MSG_CHECKING([If -DMPICH_IGNORE_CXX_SEEK is necessary])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
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
  AC_LANG_RESTORE
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
     [$2_CONFIG_OPTIONS="include $itaps_dir/lib/${itaps_pfx}-Defs.inc"],
     [AC_CHECK_FILE([$itaps_dir/${itaps_pfx}-Defs.inc],
       [$2_CONFIG_OPTIONS="include $itaps_dir/${itaps_pfx}-Defs.inc"],
       [HAVE_$2=no])])
    ;;
esac

AC_SUBST([$2_CONFIG_OPTIONS])
if test "xyes" = "x$HAVE_$2"; then 
  AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_$2"
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
AC_CHECK_FILE([$IMESH_DIR/lib/libMOAB.a], [have_moab=yes])
AM_CONDITIONAL(HAVE_MOAB, [test "xyes" = "x$have_moab"])
AC_SUBST([HAVE_MOAB])

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
MK_OPTIONAL_TOOL([algs],    [yes])
MK_OPTIONAL_TOOL([rgg],    [yes])
MK_OPTIONAL_TOOL([utils],    [yes])


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
  HAVE_CAMAL_PAVER=no
  HAVE_CAMAL_TRIADV=no
  HAVE_CAMAL_TRIDELAUNAY=no
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
    HAVE_CAMAL_SWEEP=no
    HAVE_CAMAL_TET=no
  else
    if test "x$CAMAL_VERSION" = "x510"; then
      AC_CHECK_HEADER( [CMLTriAdvance.hpp],[HAVE_CAMAL_TRIADV=yes], [HAVE_CAMAL_TRIADV=no] )
      AC_CHECK_HEADER( [CMLTriDelaunay.hpp],[HAVE_CAMAL_TRIDEL=yes], [HAVE_CAMAL_TRIDEL=no] )
      AC_CHECK_HEADER( [CMLPaver.hpp],    [HAVE_CAMAL_PAVER=yes],  [HAVE_CAMAL_PAVER=no] )
      AC_CHECK_HEADER( [CMLSweeper.hpp],  [HAVE_CAMAL_SWEEP=yes],  [HAVE_CAMAL_SWEEP=no] )
      AC_CHECK_HEADER( [CMLTetMesher.hpp],[HAVE_CAMAL_TET=yes],    [HAVE_CAMAL_TET=no] )
      AC_CHECK_HEADER( [CMLSurfEval.hpp], [HAVE_CAMAL=yes],        [HAVE_CAMAL=no] )

      CAMAL_LIBS="-lCMLSizer"
      AC_CHECK_LIB( [CMLPave],   [main], [CAMAL_LIBS="-lCMLPave $CAMAL_LIBS"], [HAVE_CAMAL_PAVER=no],  [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLTri],    [main], [CAMAL_LIBS="-lCMLTri $CAMAL_LIBS"], [HAVE_CAMAL_TRIADV=no], [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLDelTri], [main], [CAMAL_LIBS="-lCMLDelTri $CAMAL_LIBS"], [HAVE_CAMAL_TRIDEL=no], [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLSweep],  [main], [CAMAL_LIBS="-lCMLSweep $CAMAL_LIBS"], [HAVE_CAMAL_SWEEP=no],  [$CAMAL_LIBS] )
      AC_CHECK_LIB( [CMLTet],    [main], [CAMAL_LIBS="-lCMLTet $CAMAL_LIBS"], [HAVE_CAMAL_TET=no],    [$CAMAL_LIBS] )
    else
      AC_CHECK_HEADER( [CMLTriMesher.hpp],[HAVE_CAMAL_TRIADV=yes], [HAVE_CAMAL_TRIADV=no] )
      AC_CHECK_HEADER( [CMLPaver.hpp],    [HAVE_CAMAL_PAVER=yes],  [HAVE_CAMAL_PAVER=no] )
      AC_CHECK_HEADER( [CMLSweeper.hpp],  [HAVE_CAMAL_SWEEP=yes],  [HAVE_CAMAL_SWEEP=no] )
      AC_CHECK_HEADER( [CMLTetMesher.hpp],[HAVE_CAMAL_TET=yes],    [HAVE_CAMAL_TET=no] )
      AC_CHECK_HEADER( [CMLGeomEval.hpp], [HAVE_CAMAL=yes],        [HAVE_CAMAL=no] )

      CAMAL_LIBS="-lCMLShr$CAMAL_VERSION -lCMLUtl$CAMAL_VERSION -lCMLShr$CAMAL_VERSION"
      AC_CHECK_LIB( [CMLPave$CAMAL_VERSION],  [main], [CAMAL_LIBS="-lCMLPave$CAMAL_VERSION  $CAMAL_LIBS"], [HAVE_CAMAL_PAVER=no],  [$CAMAL_LIBS] )
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

AM_CONDITIONAL( HAVE_CAMAL,        [test $HAVE_CAMAL        = yes] )
AM_CONDITIONAL( HAVE_CAMAL_TRIADV, [test $HAVE_CAMAL_TRIADV = yes] )
AM_CONDITIONAL( HAVE_CAMAL_TRIDEL, [test $HAVE_CAMAL_TRIDEL = yes] )
AM_CONDITIONAL( HAVE_CAMAL_PAVER,  [test $HAVE_CAMAL_PAVER  = yes] )
AM_CONDITIONAL( HAVE_CAMAL_SWEEP,  [test $HAVE_CAMAL_SWEEP  = yes] )
AM_CONDITIONAL( HAVE_CAMAL_TET,    [test $HAVE_CAMAL_TET    = yes] )


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
#                           Output Files
################################################################################
AC_SUBST([INCLUDES])
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
                 algs/Qslim/Makefile
				 rgg/Makefile
                 utils/Makefile
				 camal/Makefile
				 plibs/Makefile
                 ])

AC_OUTPUT

AC_MSG_RESULT([C:   $CC $CFLAGS])
AC_MSG_RESULT([C++: $CXX $CXXFLAGS])
AC_MSG_RESULT([Fortran: $FC $FCFLAGS])

if test "x$IMESH_MISSING" = "xyes"; then
  AC_MSG_WARN([Could not find iMesh library.  Build will FAIL!!! Try \"--with-imesh=DIR\"])
fi
