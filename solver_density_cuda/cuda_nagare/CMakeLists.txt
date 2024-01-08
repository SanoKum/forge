find_package(CUDA REQUIRED)

if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})

    CUDA_ADD_LIBRARY(cuda_nagare STATIC
		cudaConfig.cuh
        cudaWrapper.cu
        cudaWrapper.cuh
        calcStructualVariables_d.cu
        calcStructualVariables_d.cuh
        ###boundaryCond_d.cu
        ###boundaryCond_d.cuh
        ###calcGradient_d.cu
        ###calcGradient_d.cuh
        ###calcConvDiff_d.cu
        ###calcConvDiff_d.cuh
        ###updateCenterVelocity_d.cu
        ###updateCenterVelocity_d.cuh
        ###interpVelocity_c2p_d.cu
        ###interpVelocity_c2p_d.cuh
        #solvePoisson_amgx_d.cu
        #solvePoisson_amgx_d.cuh
    )

	#target_link_libraries( amgcl_cuda cusparse )
	#link_directories(/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/11.8/targets/x86_64-linux/lib)
    #include_directories("${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
endif()
