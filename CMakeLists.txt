cmake_minimum_required (VERSION 3.5)
project (tLBM VERSION 0.1
               LANGUAGES CXX)

find_package(MPI REQUIRED)


#[==[
this is a multi-line comment put here for the hell of it.
so, here is a second line of comments...
#]==]

#include_directories ("${PROJECT_SOURCE_DIR}/include")

add_executable(tLBM src/tLBM.cpp 
                    src/Problem.cpp
                    src/TLBM_Partition.cpp
                    src/LatticeStructure.cpp
                    include/Problem.h
                    include/TLBM_Partition.h
                    include/LatticeStructure.h)

target_compile_options(tLBM PUBLIC -Wall -Wextra -Wpedantic)
target_include_directories(tLBM PUBLIC include
                                PRIVATE src ${MPI_INCLUDE_PATH})

target_link_libraries(tLBM ${MPI_C_LIBRARIES}) 
