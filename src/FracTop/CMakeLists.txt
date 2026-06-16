add_library(FracTop
    SolveRLProblem.cc
    SolveRLProblem.h
)

target_include_directories(FracTop
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(FracTop
    PRIVATE
    PUBLIC
        msg_lib
        Problem
        Operators
        utils_lib       
)
