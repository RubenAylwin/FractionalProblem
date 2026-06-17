add_library(msg_lib OBJECT
    Msg.cc
    Msg.h
)

add_library(utils_lib
    Utilities.cc
    Utilities.h
    MyTypes.cc
    Mytypes.h
)

target_link_libraries(utils_lib
    PRIVATE
        msg_lib
        Boost::boost
        Geometry
        Integration
    PUBLIC
        Eigen3::Eigen
        Discrete
)

target_include_directories(msg_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(utils_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(Discrete)
add_subdirectory(EmpiricalInterpolation)
add_subdirectory(FracTop)
add_subdirectory(Functions)
add_subdirectory(Geometry)
add_subdirectory(GreenFunc)
add_subdirectory(Integration)
add_subdirectory(Operators)
add_subdirectory(Problem)
add_subdirectory(QMC)
add_subdirectory(progressbar)
