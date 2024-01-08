cmake_minimum_required(VERSION 3.18)
project(Nagare CXX)
# hello.cppとgood_morning.cppをコンパイルして静的ライブラリlibgreetings.aを作成
# a.outという実行ファイルをmain.cppから作成

#set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE string "g++ compiler" FORCE)
set(CMAKE_CXX_COMPILER /usr/bin/g++ )

#add_definitions( -std=c++11 -Wall )
add_definitions( -std=c++17 )

set(CMAKE_CXX_FLAGS "-O3 ")

include_directories(./ /home/kumpei/src/HighFive/include /home/kumpei/util/hdf5-1.14.0/include /home/kumpei/util/amgcl/include)
link_directories(./ /home/kumpei/util/hdf5-1.14.0/lib)
add_executable(Nagare main.cpp)
target_link_libraries(Nagare hdf5 yaml-cpp)