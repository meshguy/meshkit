
#find SZIP
find_package(SZIP REQUIRED)

# add dummy target so that we can attach properties and dependencies work properly
add_dummy_external_project(szip)

#explicitly setup the include dir since add_external_project
#isn't called
set_include_dir(szip ${SZIP_INCLUDE_DIRS})

#don't use LIBRARIES or LIBRARY they are wrong
set_libraries(szip ${SZIP_LIBRARY_RELEASE})