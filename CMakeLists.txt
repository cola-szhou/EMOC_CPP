cmake_minimum_required(VERSION 3.10)
project(EMOC)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "set build type to release")    			    # make no difference in windows for Visual Studio
set(CMAKE_CONFIGURATION_TYPES "Release")	# this is for Visual Studio
set(BUILD_PYTHON_DLL OFF CACHE BOOL "Whether build python dll")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Configuration Type: ${CMAKE_CONFIGURATION_TYPES}")

# search files
file(GLOB_RECURSE SRC ./src/*.cpp)
file(GLOB_RECURSE HEADERS ./src/*.h)
file(GLOB_RECURSE PYTHON_MAIN_NAME ./src/python_main.cpp)

# set HEADERS
set(HEADERS ${HEADERS}
		vendor/random/random.h
)      

# set SRC		
set(SRC ${SRC}
		vendor/random/random.cpp
)

# set include dir and link dir
set(INCLUDE_DIR 
    ./src
    ./vendor
)

# create real EMOC project
add_subdirectory(vendor/pybind11)
pybind11_add_module(EMOC ${HEADERS} ${SRC})

# add include dir and link dir
target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})
target_link_directories(${PROJECT_NAME} PUBLIC ${LINK_DIR})