# I have no idea which versions of CMake are able to build
# the project using this specification. This is only specified
# to suppress the warning that results if it isn't.
cmake_minimum_required(VERSION 2.6)

# The project is called "mazing", and is written in C99.
#
# "-std=c99" works with gcc. If it doesn't work with your
# compiler, let me know what does and I'll try to support it.
#
# (If you're using MSVC then this is not your lucky day, because
# that compiler does not support C99 at all.)
project(mazing C)
set(CMAKE_C_FLAGS -std=c99)

# We need GMP
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_BINARY_DIR}/CMakeModules")
find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIR})

# The "exe" target builds the mazing executable
add_executable(exe main.c mazing.c fmc.c)
target_link_libraries(exe ${GMP_LIBRARIES})
set_target_properties(exe PROPERTIES OUTPUT_NAME mazing)

# Source and header files that make up the mazing library
set(LIB_SOURCES mazing.c fmc.c)
set(LIB_HEADERS mazing.h fmc.h)

# The "static" target builds the static library
add_library(static STATIC ${LIB_SOURCES})
target_link_libraries(static ${GMP_LIBRARIES})
set_target_properties(static PROPERTIES OUTPUT_NAME mazing PUBLIC_HEADER "mazing.h;fmc.h")

# The "shared" target builds the shared library
add_library(shared SHARED ${LIB_SOURCES})
target_link_libraries(shared ${GMP_LIBRARIES})
set_target_properties(shared PROPERTIES OUTPUT_NAME mazing PUBLIC_HEADER "mazing.h;fmc.h")

install(TARGETS exe RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install(TARGETS static shared
    ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mazing")
