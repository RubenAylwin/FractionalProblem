add_library(GreenFunc
    Green.cc
    GreenH.cc
    GreenL.cc
    GreenLog.cc
    GreenQP.cc
    Green.h
    GreenH.h
    GreenL.h
    GreenLog.h
    GreenQP.h
)

target_include_directories(GreenFunc
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(GreenFunc
    PRIVATE
        msg_lib
        Integration
    PUBLIC
        utils_lib       
)
