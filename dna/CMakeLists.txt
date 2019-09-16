cmake_minimum_required(VERSION 3.8)

project(cognitiv-coding-challenge-dna)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")

#set(SRC_FILES
        #person_analyzer.cpp
   #)

add_library(cogdna INTERFACE)
#add_library(cogdna SRC_FILES)
target_include_directories(cogdna
        INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_subdirectory(test)
