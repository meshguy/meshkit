add_subdirectory(core)
add_subdirectory(algs)
add_subdirectory(utils)
add_subdirectory(lemon)
add_subdirectory(extern)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp.in" "")
configure_file(
  "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp"
  @ONLY)

add_library(MeshKit
  "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp") # Force C++ language.
target_link_libraries(MeshKit
  MeshKitcore
  MeshKitalgs
  emon
  MeshKitutils
  MeshKitextern)
