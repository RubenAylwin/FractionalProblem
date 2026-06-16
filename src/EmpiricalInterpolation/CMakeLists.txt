add_library(EmpiricalInterpolation
    EmpiricalInterpolation.cc
    EmpiricalInterpolation.h
    randomized_svd.h
)

target_include_directories(EmpiricalInterpolation
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(EmpiricalInterpolation
    PRIVATE
    PUBLIC
        msg_lib
        progressbar
        utils_lib
)
