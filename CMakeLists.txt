cmake_minimum_required(VERSION 3.20)

project(VP VERSION 1.1)

file(GLOB_RECURSE SRC_FILES src/*.cc)
add_library(VP_LIB ${SRC_FILES})

find_package(Boost 1.84 COMPONENTS program_options REQUIRED)
find_package(TBB REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Charts)

target_link_libraries(VP_LIB Boost::program_options TBB::tbb Eigen3::Eigen)
target_include_directories(VP_LIB PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/includes)

target_compile_features(VP_LIB PRIVATE cxx_std_20)
target_compile_options(VP_LIB PRIVATE -Wall -Wextra -O3 -pedantic)

add_executable(VP VP.cc)
add_executable(VPQt VPQt.cc)

target_link_libraries(VP PRIVATE VP_LIB)
target_link_libraries(VPQt PRIVATE VP_LIB Qt5::Core Qt5::Widgets Qt5::Charts)