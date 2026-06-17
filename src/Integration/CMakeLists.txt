add_library(Integration
    OneDimensionalIntegration.cc
    TwoDimensionalIntegration.cc
    OneDimensionalIntegration.h
    TwoDimensionalIntegration.h
)

target_include_directories(Integration
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(Integration
    PRIVATE
        msg_lib
        Eigen3::Eigen
        Functions
    PUBLIC
        Boost::boost
        utils_lib
        
)
