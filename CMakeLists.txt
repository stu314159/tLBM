cmake_minimum_required (VERSION 3.5)
project (tLBM VERSION 0.1
               LANGUAGES CXX)

find_package(MPI REQUIRED)


#[==[
this is a multi-line comment put here for the hell of it.
so, here is a second line of comments...
#]==]

#include_directories ("${PROJECT_SOURCE_DIR}/include")

add_executable(tLBMexec src/tLBM.cpp 
                    src/Problem.cpp
                    src/TLBM_Partition.cpp
                    src/LatticeStructure.cpp
                    src/D3Q15LatticeStructure.cpp
                    src/D3Q19LatticeStructure.cpp
                    include/Problem.h
                    include/TLBM_Partition.h
                    include/LatticeStructure.h
                    include/D3Q15LatticeStructure.h
                    include/D3Q19LatticeStructure.h
                    include/TLBM_definitions.h)

target_compile_options(tLBMexec PUBLIC -Wall -Wextra -Wpedantic -std=c++11)
target_include_directories(tLBMexec PUBLIC include
	                            PRIVATE src ${MPI_CXX_INCLUDE_PATH})

target_link_libraries(tLBMexec ${MPI_CXX_LIBRARIES})

 
