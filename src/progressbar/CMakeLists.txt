add_library(progressbar INTERFACE)

target_include_directories(progressbar
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
