cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 11)

project(sav-afdist)

# Configure with cmake option: -DCMAKE_TOOLCHAIN_FILE=PATH/TO/cget.cmake
find_package(savvy REQUIRED)
find_package(Threads REQUIRED)
set(CMAKE_BUILD_TYPE Release)
add_executable(sav-afdist main.cpp)
target_link_libraries(sav-afdist savvy)
target_link_libraries(sav-afdist ${CMAKE_THREAD_LIBS_INIT})
