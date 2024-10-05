cmake_minimum_required(VERSION 3.22)
project(forge LANGUAGES CUDA)
##project(forge LANGUAGES CXX CUDA)
##
##set(CMAKE_CXX_COMPILER /usr/bin/g++ )
##add_definitions( -std=c++11 -Wall )
###add_definitions( -std=c++2a )
####set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

#cuda
set(dyn_libs_cuda cuda rt dl cudart_static cublas)

#hdf
set(hdf5_inc /home/kumpei/app/HighFive/include /home/kumpei/app/hdf/HDF5-1.14.3-Linux/HDF_Group/HDF5/1.14.3/include)
set(hdf5_libdir /home/kumpei/app/hdf/HDF5-1.14.3-Linux/HDF_Group/HDF5/1.14.3/lib)

#kdtree c++ dir
set(kdtree_dir /home/kumpei/app/kdtree/install/)

include_directories(./ ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${hdf5_inc} ${kdtree_dir}/include ) #${EIGEN_inc})
link_directories(./ ${hdf5_libdir} ${CUDA_LIB_DIR} ${kdtree_dir}/lib)

# --------------
# *** forge ***
# --------------
file( GLOB forge_cppfiles 
    main.cpp ,
    gradient.cpp ,
    mesh/mesh.cpp ,
    output/output.cpp ,
    update.cpp ,
    variables.cpp ,
    boundaryCond.cpp ,
    implicitCorrection.cpp ,
    convectiveFlux.cpp ,
    input/solverConfig.cpp ,
    input/calcWallDistance_kdtree.cpp ,
    setDT.cpp ,
    setStructualVariables.cpp ,
    dependentVariables.cpp 
)

set_source_files_properties(${forge_cppfiles} PROPERTIES LANGUAGE CUDA)

add_executable(forge ${forge_cppfiles})


add_subdirectory(cuda_forge )


target_compile_options(forge PUBLIC
             $<$<COMPILE_LANGUAGE:CXX>:-O3 -fopenmp>)


target_link_libraries(forge PRIVATE cuda_forge)
target_link_libraries(forge PRIVATE hdf5 yaml-cpp ${dyn_libs_cuda} kdtree)

# ----------------------
# *** mesh conversion***
# ----------------------
file( GLOB mesh_cppfiles 
    mesh/convertGmshToForge.cpp 
    mesh/mesh.cpp 
    variables.cpp 
    boundaryCond.cpp 
    input/solverConfig.cpp 
    input/calcWallDistance_kdtree.cpp 
)

set_source_files_properties(${mesh_cppfiles} PROPERTIES LANGUAGE CUDA)

add_executable(convertGmshToForge  ${mesh_cppfiles})
target_link_libraries(convertGmshToForge hdf5 yaml-cpp)
target_link_libraries(convertGmshToForge cuda_forge kdtree )
