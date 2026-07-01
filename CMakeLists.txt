
cmake_minimum_required(VERSION 3.20)

project(FracProb
        VERSION 2.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost 1.74 COMPONENTS program_options REQUIRED)
find_package(TBB REQUIRED)
find_package(Eigen3 REQUIRED)

add_subdirectory(src)

add_executable(VP VP.cc)
target_link_libraries(VP
    PRIVATE
        Discrete
        FracTop
        Functions
        Operators
        Problem
        progressbar
        Boost::program_options
        utils_lib
        msg_lib
)

option(TESTING "Building tests" ON)

if (TESTING)
   enable_testing()
   add_subdirectory(test)
endif()
