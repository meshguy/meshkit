find_package(Freetype)

find_path(
        FTGL_INCLUDE_DIRS
        FTGL/ftgl.h
        HINTS ${FTGL_HOME}/include
        )

find_library(
        FTGL_LIBRARIES
        ftgl
        HINTS ${FTGL_HOME}/lib
        )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
        FTGL
        DEFAULT_MSG
        FTGL_INCLUDE_DIRS
        FTGL_LIBRARIES
        )