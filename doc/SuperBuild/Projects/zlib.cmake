# zlib supports cmake. the only problem is that we need to remove the zconf.h
# file.
add_external_project(
  zlib
  # remove the zconf.h as a patch step.
  PATCH_COMMAND ${CMAKE_COMMAND} -E remove -f <SOURCE_DIR>/zconf.h
  )
ExternalProject_Get_Property(${name} install_dir)
set_libraries(zlib ${install_dir}/lib/libz${CMAKE_SHARED_LIBRARY_SUFFIX})
