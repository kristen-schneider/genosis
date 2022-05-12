cmake_minimum_required(VERSION 3.22)
project("precision-medicine" CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-local")
include(ExternalProject)
set(CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "-std=c++0x")


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

include_directories ("${PROJECT_SOURCE_DIR}/include" )

include_directories(src)

add_executable (precision-medicine "src/main.cpp")
