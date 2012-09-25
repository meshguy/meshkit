
#find ftgl
find_package(FTGL REQUIRED)

# add dummy target so that we can attach properties and dependencies work properly
add_dummy_external_project(ftgl)

#explicitly setup the include dir since add_external_project
#isn't called
set_include_dir(ftgl ${FTGL_INCLUDE_DIRS})
