#Open Cascade Community Edition supports CMake
add_external_project(OCE
  DEPENDS ftgl  
  CMAKE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=TRUE
    -OCE_DISABLE_X11:BOOL=FALSE
    -OCE_VISU_DEP:BOOL=TRUE
    -OCE_VISUALISATION:BOOL=TRUE
    -OCE_DRAW:BOOL=TRUE
    -DFTGL_INCLUDE_DIR:PATH=<INSTALL_DIR>/FTGL
    -DOCE_INSTALL_PREFIX:FilePath=<INSTALL_DIR>
    #force the include dir path, so it doesn't install the include files into
    #install/include/oce/, because than CGM can't find it.
    -DOCE_INSTALL_INCLUDE_DIR:FilePath=include
)

# remove the installed oce-config.h from the install tree
# so that the build doesn't use that over the configured header in the build
# directory. If we used the installed oce-config the incremental builds will
# fail as it doesn't have all the defines as the build version
ExternalProject_Add_Step(OCE OCE-remove-config-file
    COMMAND  ${CMAKE_COMMAND} -E remove -f <INSTALL_DIR>/include/oce-config.h
    COMMENT "Removing the installed oce-config.h so that incremental builds works."
    DEPENDEES configure
    DEPENDERS build
    ALWAYS 1
    )