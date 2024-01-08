cmake_minimum_required(VERSION 3.18)
project(nagare CXX)
# hello.cppとgood_morning.cppをコンパイルして静的ライブラリlibgreetings.aを作成
# a.outという実行ファイルをmain.cppから作成

#set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE string "g++ compiler" FORCE)
set(CMAKE_CXX_COMPILER /usr/bin/g++ )

#add_definitions( -std=c++11 -Wall )
add_definitions( -std=c++2a )

set(CMAKE_CXX_FLAGS "-O3")

include_directories(./ /home/kumpei/src/HighFive/include /home/kumpei/util/hdf5-1.14.0/include /home/kumpei/util/amgcl/include)
link_directories(./ /home/kumpei/util/hdf5-1.14.0/lib)

#add_library(nagare STATIC gradient.cpp)

add_executable(nagare main.cpp gradient.cpp mesh/mesh.cpp output/output.cpp update.cpp variables.cpp boundaryCond.cpp solvePoisson_amgcl.cpp solveNavierStokes.cpp input/solverConfig.cpp calcCFL.cpp setStructualVariables.cpp)

add_executable(convertGmshToNagare mesh/convertGmshToNagare.cpp mesh/mesh.cpp variables.cpp boundaryCond.cpp input/solverConfig.cpp )

target_link_libraries(nagare hdf5 yaml-cpp)
target_link_libraries(convertGmshToNagare hdf5 yaml-cpp)

add_subdirectory(amgcl_cuda)
target_include_directories(nagare PUBLIC amgcl_cuda)
target_link_libraries(nagare amgcl_cuda)

