cmake_minimum_required(VERSION 2.8.3)
project(graph_slam)
# set(CMAKE_BUILD_TYPE Release)
## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++14)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED
  roscpp
  rospy
  std_msgs
  tf
  sensor_msgs
  geometry_msgs
  laser_geometry
)

find_package(PCL 1.7 REQUIRED)

list( APPEND CMAKE_MODULE_PATH /home/touchair/ekf_ws/src/slambook2/3rdparty/g2o/cmake_modules ) 
set(G2O_ROOT /usr/local/include/g2o) 

find_package(G2O REQUIRED)



catkin_package(
 INCLUDE_DIRS include
#  LIBRARIES gaussian_newton_scanmatcher
#  CATKIN_DEPENDS other_catkin_pkg
#  DEPENDS system_lib
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
 include
 ${EIGEN3_INCLUDE_DIR}
 ${catkin_INCLUDE_DIRS}
 ${PCL_INCLUDE_DIRS}
 ${G2O_INCLUDE_DIR}
 ${G2O_INCLUDE_DIRS}
)

link_directories(${G2O_LIBRARY_DIRS})

## Declare a C++ library
# add_library(${PROJECT_NAME}
#   src/${PROJECT_NAME}/gaussian_newton_scanmatcher.cpp
# )

add_library(${PROJECT_NAME}
  src/gaussian_newton_method.cpp
  src/keyframe.cpp
  src/map_cspace.cpp
  src/map.c
)

## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Declare a C++ executable
## With catkin_make all packages are built within a single CMake context
## The recommended prefix ensures that target names across packages don't collide
# add_executable(gaussian_newton_node src/gaussian_newton_method.cpp src/map.c src/map_cspace.cpp src/main.cpp src/graph_slam.cpp src/back_end.cpp src/keyframe.cpp)

add_executable(front_end src/front_end.cpp)

## Specify libraries to link a library or executable target against
target_link_libraries(front_end
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
  ${PCL_LIBRARIES}
  ${G2O_CORE_LIBRARY}
  ${G2O_STUFF_LIBRARY}
)

add_executable(back_end src/back_end.cpp)

target_link_libraries(back_end
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
  ${PCL_LIBRARIES}
  ${G2O_CORE_LIBRARY}
  ${G2O_STUFF_LIBRARY}
)



#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_gaussian_newton_scanmatcher.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
