cmake_minimum_required(VERSION 3.12)
project(fluidanimate)

set(CMAKE_CXX_STANDARD 23)
# set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3")

add_executable(${PROJECT_NAME} pthreads.cpp cellpool.hpp cellpool.cpp parsec_barrier.hpp 
parsec_barrier.cpp cha.hpp cha.cpp tile.hpp topology.hpp topology.cpp)

target_link_libraries(${PROJECT_NAME} pthread)
target_compile_definitions(${PROJECT_NAME} PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)