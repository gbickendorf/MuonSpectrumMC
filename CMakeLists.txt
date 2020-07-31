# only for cmake --version >= 3.5.1
cmake_minimum_required(VERSION 3.5.1)

# project name
project(test)

# creates the variable EXEC and sets it to hello
set(EXEC BSM-MuonSpectrum)

# set the C++14 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# I../includes
include_directories(includes)


# puts all .cpp files inside src to the SOURCES variable
file(GLOB SOURCES src/*.cpp)

# compiles the files defined by SOURCES to generante the executable defined by EXEC
add_executable(${EXEC} ${SOURCES})

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(${EXEC} GSL::gsl GSL::gslcblas)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${EXEC} OpenMP::OpenMP_CXX)
endif()
