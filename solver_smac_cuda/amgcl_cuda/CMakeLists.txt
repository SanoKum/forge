find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})
    CUDA_ADD_LIBRARY(amgcl_cuda STATIC
        solvePoisson_amgcl_cuda.cu
        solvePoisson_amgcl_cuda.h
    )
	target_link_libraries( amgcl_cuda ${CUDA_cusparse_LIBRARY} )
endif()
