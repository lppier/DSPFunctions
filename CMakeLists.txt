cmake_minimum_required(VERSION 3.4)
project(DSPFunctions)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories("/usr/include/eigen3/Eigen")
include_directories("/usr/local/include")

link_directories("/opt/intel/compilers_and_libraries_2016.0.109/linux/ipp/lib/intel64")
link_directories("/usr/local/lib")

link_libraries(ipps)
link_libraries(ippi)
link_libraries(ippcore)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "/home/pier/gitworkspace/DSPFunctions/Debug")

set(SOURCE_FILES
        DSPFunctions.hpp
        DSPFunctions.cpp
        MappingFunctions.cpp
        MappingFunctions.hpp
        MovingAvgFilter.h
        main.cpp)
add_executable(DSPFunctions ${SOURCE_FILES})