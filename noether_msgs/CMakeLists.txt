cmake_minimum_required(VERSION 3.5)
project(noether_msgs)

set(MESSAGE_DEPENDENCIES geometry_msgs)
set(DEPENDENCIES ament_cmake rosidl_default_generators ${MESSAGE_DEPENDENCIES})
foreach(DEP ${DEPENDENCIES})
  find_package(${DEP} REQUIRED)
endforeach()

rosidl_generate_interfaces(${PROJECT_NAME}
  "msg/ToolPath.msg"
  "msg/ToolPaths.msg"
  DEPENDENCIES ${MESSAGE_DEPENDENCIES})

ament_export_dependencies(${DEPENDENCIES})

ament_package()
