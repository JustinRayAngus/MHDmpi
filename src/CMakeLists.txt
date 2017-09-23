project("MHDwithMPI")
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# !!! PHYScpp and PHYS_HOME must be set as environment 
# variables prior to "cmake pathtoCMakeLists.txt"
# i.e
# export PHYScpp=burgers1D.cpp
# export PHYS_HOME="../burgers1D"
#set(PHYScpp burgers1D.cpp)
#set(ENV{PHYS_HOME} "../burgers1D")
include_directories($ENV{PHYS_HOME})

# include stuff for HDF5
set(ENV{HDF5_HOME} "/Users/angus1/HDF5/hdf5-1.10.0-patch1/hdf5")
#find_package(HDF5 "1.8.14" REQUIRED COMPONENTS C HL)
#include_directories(${HDF5_INCLUDE_DIRS})
#link_directories(${HDF5_LIBRARY_DIRS})
include_directories($ENV{HDF5_HOME}/include)
link_directories($ENV{HDF5_HOME}/lib)

# include stuff for mpi
#set(ENV{MPI_HOME} "/opt/local")
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
#include_directories($ENV{MPI_HOME}/include)
#link_directories($ENV{MPI_HOME}/lib)

# include path to header files
# (needed for #include in PHYScpp)
#set(ENV{MHDmpi_HOME} "/Users/angus1/Programs/MHDmpi")
#include_directories($ENV{MHDmpi_HOME}/src)
include_directories(./)

add_definitions(-std=c++0x)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

set(sources 
   main.cpp
   domainGrid.cpp
   vectorMath.cpp
   timeDomain.cpp
   HDF5dataFile.cpp
   jsoncpp.cpp
   $ENV{PHYScpp})
#set(sources ${sources} ${PHYScpp})

add_executable(mhd ${sources})
target_link_libraries(mhd ${MPI_LIBRARIES} hdf5 hdf5_cpp)


