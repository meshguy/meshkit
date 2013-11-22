
#find zlib
find_package(ZLIB REQUIRED)

# add dummy target so that we can attach properties and dependencies work properly
add_dummy_external_project(zlib)

#explicitly setup the include dir since add_external_project
#isn't called
set_include_dir(zlib ${ZLIB_INCLUDE_DIRS})
set_libraries(zlib ${ZLIB_LIBRARIES})