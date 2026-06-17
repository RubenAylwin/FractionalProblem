add_executable(test_fractional
    test.cc
)

target_include_directories(test_fractional
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_fractional
    PRIVATE
        msg_lib
        utils_lib
        Discrete
        EmpiricalInterpolation
        FracTop
        Functions
        Geometry
        GreenFunc
        Integration
        Operators
        Problem
        QMC
        progressbar
)


add_test(NAME Tests COMMAND test_fractional)
