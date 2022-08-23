cmake_minimum_required(VERSION 3.5)
project(MSbrkga)

include_directories(/opt/boost_1_79_0)

link_directories(/opt/boost_1_79_0)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
	BRKGA.h
	MSDecoder.h
	Population.h
	MTRand.h
	samplecode.cpp
	MSDecoder.cpp)

add_executable(MSbrkga ${SOURCE_FILES})
