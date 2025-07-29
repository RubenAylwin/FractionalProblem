cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra  -Xpreprocessor -O3 -pedantic")

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

project(hallo VERSION 1.0)
file(GLOB_RECURSE SRC_FILES src/*.cc)

find_package(Boost 1.84 COMPONENTS program_options REQUIRED)
find_package(TBB REQUIRED)
add_executable(VP VP.cc ${SRC_FILES})
target_include_directories(VP PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/includes ~/eigen)
target_link_libraries(VP Boost::program_options TBB::tbb)
