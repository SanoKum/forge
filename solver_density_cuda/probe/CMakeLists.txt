project(probe LANGUAGES CUDA)

add_library(probe STATIC
    point_probes.cuh
    point_probes.cu
)

target_include_directories(probe PRIVATE .. )

set_property(TARGET probe  PROPERTY CUDA_ARCHITECTURES 86)
target_compile_options(probe PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_86,code=sm_86> -lineinfo)


