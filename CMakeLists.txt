cmake_minimum_required (VERSION 2.8)

project (MDAtom)

include_directories(
	include
)

SET (HEADERS
	include/mdprog.h
)
 
SET (SOURCE
	src/main.cpp
	src/runmda.cpp
	src/auxiliary.cpp
	src/cenma.cpp
	src/confa.cpp
	src/forcea.cpp
	src/gauss.cpp
	src/packa.cpp
	src/random.cpp
	src/shia.cpp
)

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

add_executable(mdatom ${SOURCE} ${HEADERS})
