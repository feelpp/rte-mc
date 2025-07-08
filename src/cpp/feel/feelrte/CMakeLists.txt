add_library(feelrte
  SHARED
  geometry.cpp
)
target_include_directories(feelrte
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(feelrte
  PUBLIC
    Feelpp::feelpp       # Feel++ core
)

