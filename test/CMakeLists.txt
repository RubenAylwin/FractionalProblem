cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-O3")

project(hallo VERSION 1.0)
file(GLOB_RECURSE SRC_FILES ../src/*.cc)

find_package(Boost 1.84 COMPONENTS program_options unit_test_framework REQUIRED)
find_package(TBB REQUIRED)
find_package(Eigen3 REQUIRED)
add_executable(test test.cc ${SRC_FILES})
target_include_directories(test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../includes)
target_link_libraries(test PRIVATE Boost::program_options TBB::tbb Eigen3::Eigen)
target_compile_options(test PUBLIC -Wall -Xpreprocessor -pedantic)
