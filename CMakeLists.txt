cmake_minimum_required(VERSION 3.5)
project(MSbrkga)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++0x ")

include_directories(/opt/boost_1_73_0)

link_directories(/opt/boost_1_73_0)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
	BRKGA.h
	MSDecoder.h
	Population.h
	MTRand.h
	SampleDecoder.h
	SampleDecoder.cpp
	samplecode.cpp
	MSDecoder.cpp main.cpp)

add_executable(MSbrkga ${SOURCE_FILES})
