cmake_minimum_required(VERSION 3.2.2)
project(gp_comet)


set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

add_subdirectory(external_library)
add_subdirectory(library)
add_subdirectory(app)
