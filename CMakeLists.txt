cmake_minimum_required(VERSION 2.8.4)
project(pcluster)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    main.cc
    main.h
    nw.cc
    nw.h
    nw_result.cc)

add_executable(pcluster ${SOURCE_FILES})