
#find OCE
find_package(OCE REQUIRED)

# add dummy target so that we can attach properties and dependencies work properly
add_dummy_external_project(OCE)

#explicitly setup the include dir since add_external_project
#isn't called
set_include_dir(OCE ${OCE_INCLUDE_DIRS})
