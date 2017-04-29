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

#add_library(MeshKit
#  "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp") # Force C++ language.

link_libraries(MeshKit
#  PRIVATE
    MeshKitcore
    MeshKitalgs
    lemon
    MeshKitutils
    MeshKitextern)

add_library(MeshKit
     $<TARGET_OBJECTS:MeshKitcore>
      $<TARGET_OBJECTS:MeshKitalgs>
      $<TARGET_OBJECTS:MeshKitutils>
      $<TARGET_OBJECTS:lemon>
      $<TARGET_OBJECTS:MeshKitextern>
      )

install(
  TARGETS  MeshKit
  EXPORT   MeshKitTargets
  RUNTIME  DESTINATION bin
  LIBRARY  DESTINATION lib
  ARCHIVE  DESTINATION lib
  COMPONENT runtime)

set_property(GLOBAL APPEND
  PROPERTY meshkit_export_targets
  MeshKit)
