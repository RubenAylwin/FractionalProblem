add_library(Functions
    ParametrizedFunction.cc
    ScalarValuedFunction.cc
    ParametrizedFunction.h
    ScalarValuedFunction.h
)

target_include_directories(Functions
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(Functions
    PRIVATE
        msg_lib
    PUBLIC
        utils_lib
)
