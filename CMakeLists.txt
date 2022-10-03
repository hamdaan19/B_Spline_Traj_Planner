cmake_minimum_required(VERSION 3.0.2)
project(B-Spline-Traj-Planner)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
  mavros
  roscpp
  rospy
  std_msgs
  geometry_msgs
  trajectory_msgs
  mavros_msgs
  message_generation
  ompl
)

## System dependencies are found with CMake's conventions
find_package(ompl REQUIRED)
find_package(Eigen3 3.1.0 REQUIRED NO_MODULE)
find_package(nlopt REQUIRED)
find_package(GSL REQUIRED)

add_service_files(
  FILES
  TrajPlanner.srv
)

generate_messages(
  DEPENDENCIES
  geometry_msgs
  trajectory_msgs
  mavros_msgs
  std_msgs
)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES mav_companion_api
   CATKIN_DEPENDS mavros roscpp rospy std_msgs message_runtime
#  DEPENDS system_lib
)

set(PYTHON_INCLUDE_DIR "/usr/include/python3.8")

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${OMPL_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include
  ${Boost_INCLUDE_DIR}
  ${EIGEN3_INCLUDE_DIR}
  ${NLOPT_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIR}
)

add_executable(
  trajPlannerServer
  src/traj_planner.cpp
  src/traj_planner_server.cpp
  src/BSpline.cpp
  src/TrajOpt.cpp
)

target_link_libraries(
    trajPlannerServer
    ${OMPL_LIBRARIES}
    ${catkin_LIBRARIES}
    ${EIGEN3_LIBS}
    ${NLOPT_LIBRARIES}
    -lnlopt
    -lm
    -lgsl 
    -lgslcblas
    -lpython3.8
)

add_dependencies(trajPlannerServer B-Spline-Traj-Planner_gencpp)