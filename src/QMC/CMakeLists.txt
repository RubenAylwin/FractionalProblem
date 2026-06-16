add_library(QMC
    Halton.cc
    Sobol.cc
    Halton.h
    Sobol.h
)

target_include_directories(QMC
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
