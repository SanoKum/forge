project(probe LANGUAGES CUDA)

add_library(probe STATIC
    point_probes.cuh
    point_probes.cu
)

target_include_directories(probe PRIVATE .. )

if(CMAKE_CUDA_ARCHITECTURES)
    set_property(TARGET probe PROPERTY CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
endif()
target_compile_options(probe PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3 -lineinfo>)


