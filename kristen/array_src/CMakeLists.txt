cmake_minimum_required(VERSION 3.21)
project(precision-medicine)

set(CMAKE_CXX_STANDARD 14)

include_directories("${PROJECT_SOURCE_DIR}")
include_directories(
		${PROJECT_SOURCE_DIR}/include
		${PROJECT_SOURCE_DIR}/src
)

file(GLOB all_SRCS
        "${PROJECT_SOURCE_DIR}/src/main.cpp"
        )

add_executable(precision-medicine ${all_SRCS})
