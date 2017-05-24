add_subdirectory(extern)
add_subdirectory(utils)
add_subdirectory(lemon)
add_subdirectory(algs)
add_subdirectory(core)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp.in" "")
configure_file(
  "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp"
  @ONLY)

#add_library(MeshKit
#  "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp") # Force C++ language.

add_library(MeshKit
    ${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp
#     MeshKitcore
#      MeshKitalgs
#      MeshKitutils
#      lemon
#      MeshKitextern
)

target_link_libraries(MeshKit
  PUBLIC
  MeshKitAlgs   
  MeshKitextern    
  MeshKitutils        
  MeshKitcore    
  lemon
)


install(
  TARGETS  MeshKit
  EXPORT   MeshKit
  RUNTIME  DESTINATION bin
  LIBRARY  DESTINATION lib
  ARCHIVE  DESTINATION lib
  COMPONENT runtime)

#set_property(GLOBAL APPEND
#  PROPERTY meshkit_export_targets
#  MeshKit)

get_property(meshkit_targets GLOBAL
  PROPERTY meshkit_export_targets)

export(
      TARGETS ${meshkit_targets}
      FILE    "${CMAKE_BINARY_DIR}/MeshKitTargets.cmake")
    
    if (WIN32)
      set(cmakedir cmake)
    else ()
      set(cmakedir lib/cmake/MeshKit)
    endif ()


install(
    EXPORT      MeshKit
    DESTINATION "${cmakedir}"
    FILE        MeshKitTargets.cmake
    COMPONENT   development)
  #EXPORT      MESHKITTargets
  #DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOAB/)

configure_file(
  "${CMAKE_SOURCE_DIR}/MeshKitConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/MeshKitConfig.cmake"
  @ONLY)
install(
  FILES       "${CMAKE_BINARY_DIR}/MeshKitConfig.cmake"
  DESTINATION "${cmakedir}"
  COMPONENT   development)
