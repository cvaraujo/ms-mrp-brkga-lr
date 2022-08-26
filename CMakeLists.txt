cmake_minimum_required(VERSION 3.5)
project(MSbrkga)

find_package(OpenMP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++0x ")

include_directories(/opt/boost_1_79_0)

link_directories(/opt/boost_1_79_0)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
	BRKGA.h
	MSDecoder.h
	Population.h
	MTRand.h
	samplecode.cpp
	MSDecoder.cpp
	edmonds_optimum_branching.hpp
    edmonds_optimum_branching_impl.hpp)

add_executable(MSbrkga ${SOURCE_FILES} $<$<CONFIG:Release>:OpenMP::OpenMP_CXX>)
