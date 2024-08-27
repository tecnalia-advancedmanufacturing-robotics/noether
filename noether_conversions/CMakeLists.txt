cmake_minimum_required(VERSION 3.5)
project(noether_conversions)

# Add support for C++11
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()

set(AMENT_DEPENDENCIES
  noether_msgs
  geometry_msgs
  visualization_msgs
  shape_msgs
  tf2_eigen
)
set(DEPENDENCIES
  ament_cmake ${AMENT_DEPENDENCIES}
)

foreach(DEP ${DEPENDENCIES})
  find_package(${DEP} REQUIRED)
endforeach()

set(CMAKE_CXX_STANDARD 14)

find_package(VTK REQUIRED NO_MODULE)
if(VTK_FOUND AND "${VTK_VERSION}" VERSION_LESS 7.1)
  message(FATAL_ERROR
    "The minimum required version of VTK is 7.1, but found ${VTK_VERSION}")
  set(VTK_FOUND FALSE)
else()
  include(${VTK_USE_FILE})
endif()

find_package(PCL REQUIRED COMPONENTS IO)
if(PCL_FOUND AND "${PCL_VERSION}" VERSION_LESS 1.9)
  message(WARNING
    "The minimum required version of PCL is 1.9, but found ${PCL_VERSION} in path first. Checking for exactly 1.9")
  find_package(PCL REQUIRED 1.9)
else()
  add_definitions(${PCL_DEFINITIONS})
endif()

find_package(Eigen3 REQUIRED)

set(INCLUDE_DIRS include ${ament_cmake_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS} ${PCL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
include_directories(${INCLUDE_DIRS})


add_library(${PROJECT_NAME} src/${PROJECT_NAME}.cpp)
target_link_libraries(${PROJECT_NAME} ${LIBS} ${VTK_LIBRARIES} ${PCL_LIBRARIES})
ament_target_dependencies(${PROJECT_NAME} ${AMENT_DEPENDENCIES})

#############
## Install ##
#############
install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib/${PROJECT_NAME})

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION include/${PROJECT_NAME}
  PATTERN ".svn"
  EXCLUDE)

ament_export_dependencies(${DEPENDENCIES} VTK PCL)
ament_export_include_directories(${INCLUDE_DIRS})
ament_export_libraries(${PROJECT_NAME} ${LIBS})

ament_package()
