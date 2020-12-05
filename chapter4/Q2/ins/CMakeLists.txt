project(ins)
cmake_minimum_required(VERSION 3.10)
add_subdirectory(${PROJECT_SOURCE_DIR}/GeographicLib)
include_directories(${PROJECT_SOURCE_DIR}/GeographicLib/include/)
include_directories(/usr/include/eigen3)

find_package(PkgConfig REQUIRED)
pkg_check_modules(YAML_CPP REQUIRED yaml-cpp)
include_directories(${YAML_CPP_INCLUDEDIR})

add_executable(ins ins.cpp)
target_link_libraries( ins libGeographiccc ${YAML_CPP_LIBRARIES})

add_executable(eskf_ins eskf_ins.cpp)
target_link_libraries( eskf_ins libGeographiccc ${YAML_CPP_LIBRARIES})
