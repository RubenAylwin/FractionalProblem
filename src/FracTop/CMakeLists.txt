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
        Problem
        msg_lib
        Operators
    PUBLIC
        utils_lib       
)
