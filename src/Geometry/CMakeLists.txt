add_library(Geometry
    Curve.cc
    GeometricVector.cc
    Point.cc
    Curve.h
    GeometricVector.h
    Point.h
)

target_include_directories(Geometry
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(Geometry
    PRIVATE
        msg_lib
        Functions
    PUBLIC
        utils_lib
)

