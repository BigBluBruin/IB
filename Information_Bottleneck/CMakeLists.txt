cmake_minimum_required(VERSION 3.1.0)
project(noise_testing) 

#Bring the headers, such as Student.h into the project
include_directories(include) 

#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp) 


#set C++ 11 standard
set(CMAKE_CXX_STANDARD 11)

#set compile flags
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE release)
endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp") 


add_executable(runner ${SOURCES})