cmake_minimum_required(VERSION 2.6)

project(mazing C)
set(CMAKE_C_FLAGS -std=c99)

find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARIES NAMES gmp libgmp)
include_directories(${GMP_INCLUDE_DIR})

add_executable(mazing mazing.c main.c fmc.c)
target_link_libraries(mazing ${GMP_LIBRARIES})
