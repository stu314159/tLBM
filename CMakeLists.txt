cmake_minimum_required (VERSION 3.17)
project (tLBM VERSION 0.2
	LANGUAGES CUDA CXX)

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(CUDA REQUIRED)


#[==[
this is a multi-line comment put here for the hell of it.
so, here is a second line of comments...
#]==]

#include_directories ("${PROJECT_SOURCE_DIR}/include")

add_executable(tLBMexec src/tLBM.cpp 
                    src/Problem.cpp
                    src/TLBM_Partition.cpp
                    include/Problem.h
                    include/TLBM_Partition.h
                    include/LatticeStructure.hpp
                    include/D3Q15LatticeStructure.hpp
                    include/D3Q19LatticeStructure.hpp
                    include/D3Q27LatticeStructure.hpp
                    include/HaloDataObject.hpp
                    include/HaloDataOrganizer.hpp
                    include/TLBM_definitions.h)

target_compile_options(tLBMexec PUBLIC -O3 -Wall -Wextra -Wpedantic -std=c++11)
target_include_directories(tLBMexec PUBLIC include
	                            PRIVATE src ${MPI_CXX_INCLUDE_PATH} ${HDF5_INCLUDE_DIRS})

target_link_libraries(tLBMexec ${MPI_CXX_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

 
