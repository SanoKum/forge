project(cuda_forge LANGUAGES CUDA)

add_library(cuda_forge STATIC
    cudaConfig.cuh
    cudaWrapper.cu
    cudaWrapper.cuh
    calcStructualVariables_d.cu
    calcStructualVariables_d.cuh
    dependentVariables_d.cu
    dependentVariables_d.cuh
    update_d.cu
    update_d.cuh
    boundaryCond_d.cu
    boundaryCond_d.cuh
    calcGradient_d.cu
    calcGradient_d.cuh
    convectiveFlux_d.cu
    convectiveFlux_d.cuh
    ducrosSensor_d.cu
    ducrosSensor_d.cuh
 
    #AUSM_d.cu
    #AUSM_d.cuh
    viscousFlux_d.cu
    viscousFlux_d.cuh

    implicitCorrection_d.cu
    implicitCorrection_d.cuh
    timeIntegration_d.cu
    timeIntegration_d.cuh
    setDT_d.cu
    setDT_d.cuh
    #derived_atomic_functions.cu
    #derived_atomic_functions.cuh

)


target_include_directories(cuda_forge PRIVATE .. )

set_property(TARGET cuda_forge PROPERTY CUDA_ARCHITECTURES 86)
target_compile_options(cuda_forge PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_86,code=sm_86>)

#set_target_properties(cuda_forge PROPERTIES CUDA_ARCHITECTURES "86")
#target_compile_options(cuda_forge PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_86,code=sm_86> -dc -w)
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -std=c++11 -gencode arch=compute_86,code=sm_86)

#target_link_libraries( amgcl_cuda cusparse )
#link_directories(/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/11.8/targets/x86_64-linux/lib)
