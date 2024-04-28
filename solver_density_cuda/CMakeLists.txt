cmake_minimum_required(VERSION 3.22)
project(forge LANGUAGES CXX CUDA)

set(CMAKE_CXX_COMPILER /usr/bin/g++ )
add_definitions( -std=c++11 -Wall )
#add_definitions( -std=c++2a )
##set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

#cuda
set(dyn_libs_cuda cuda rt dl cudart_static cublas)

#hdf
set(hdf5_inc /home/kumpei/app/HighFive/include /home/kumpei/app/hdf/HDF5-1.14.3-Linux/HDF_Group/HDF5/1.14.3/include)
set(hdf5_libdir /home/kumpei/app/hdf/HDF5-1.14.3-Linux/HDF_Group/HDF5/1.14.3/lib)

include_directories(./ ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${hdf5_inc} ) #${EIGEN_inc})
link_directories(./ ${hdf5_libdir} ${CUDA_LIB_DIR})

# --------------
# *** forge ***
# --------------
add_executable(forge main.cpp gradient.cpp mesh/mesh.cpp output/output.cpp update.cpp 
               variables.cpp boundaryCond.cpp implicitCorrection.cpp
               convectiveFlux.cpp input/solverConfig.cpp 
               setDT.cpp setStructualVariables.cpp  #solvePoisson_amgx.cpp 
               dependentVariables.cpp)

add_subdirectory(cuda_forge )


target_compile_options(forge PUBLIC
             $<$<COMPILE_LANGUAGE:CXX>:-O3 -fopenmp>)


target_link_libraries(forge PRIVATE cuda_forge)
target_link_libraries(forge PRIVATE hdf5 yaml-cpp ${dyn_libs_cuda})

# ----------------------
# *** mesh conversion***
# ----------------------
add_executable(convertGmshToForge  mesh/convertGmshToForge.cpp mesh/mesh.cpp variables.cpp boundaryCond.cpp input/solverConfig.cpp )
target_link_libraries(convertGmshToForge hdf5 yaml-cpp)
target_link_libraries(convertGmshToForge cuda_forge )
