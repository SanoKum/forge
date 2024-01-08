cmake_minimum_required(VERSION 3.18)
project(nagare CUDA CXX)


set(CMAKE_CXX_COMPILER /usr/bin/g++ )
#add_definitions( -std=c++11 -Wall )
add_definitions( -std=c++2a )
set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

#cuda
set(libs_all_cuda cusparse cusolver rt dl)
#set(dyn_libs_cuda amgxsh rt dl cudart_static cublas)
set(dyn_libs_cuda rt dl cudart_static cublas)

#amgx
#set(AMGX_INCLUDES ${THRUST_DIR} /home/kumpei/app/AMGX/include )
#set(AMGX_LIB_DIR /home/kumpei/app/AMGX/build)

#eigen
set(EIGEN_inc /usr/include/eigen3 )

#hdf
set(hdf5_inc /home/kumpei/app/HighFive/include /home/kumpei/app/hdf/HDF5-1.14.3-Linux/HDF_Group/HDF5/1.14.3/include)
set(hdf5_libdir /home/kumpei/app/hdf/HDF5-1.14.3-Linux/HDF_Group/HDF5/1.14.3/lib)

#include_directories(./ ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${hdf5_inc} ${AMGX_INCLUDES} ${EIGEN_inc})
include_directories(./ ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${hdf5_inc} ${EIGEN_inc})
#link_directories(./ ${hdf5_libdir} ${AMGX_LIB_DIR} ${CUDA_LIB_DIR})
link_directories(./ ${hdf5_libdir} ${CUDA_LIB_DIR})

# --------------
# *** nagare ***
# --------------
add_executable(nagare main.cpp gradient.cpp mesh/mesh.cpp output/output.cpp update.cpp 
               variables.cpp boundaryCond.cpp timeIntegration.cpp implicitCorrection.cpp
               convectiveFlux.cpp input/solverConfig.cpp 
               setDT.cpp setStructualVariables.cpp  #solvePoisson_amgx.cpp 
               dependentVariables.cpp)

target_include_directories(nagare PUBLIC cuda_nagare )

add_subdirectory(cuda_nagare)
target_link_libraries(nagare cuda_nagare)
target_link_libraries(nagare hdf5 yaml-cpp ${dyn_libs_cuda})
#target_link_libraries(nagare ${libs_all_cuda})

# ----------------------
# *** mesh conversion***
# ----------------------
add_executable(convertGmshToNagare mesh/convertGmshToNagare.cpp mesh/mesh.cpp variables.cpp boundaryCond.cpp input/solverConfig.cpp )
target_link_libraries(convertGmshToNagare hdf5 yaml-cpp)
target_link_libraries(convertGmshToNagare cuda_nagare)
