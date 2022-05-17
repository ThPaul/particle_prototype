cmake_minimum_required(VERSION 3.10)
project(particle_timing)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++17")
find_package(HPX REQUIRED)
add_executable(single_cell single_cell.cpp)
target_link_libraries(single_cell HPX::hpx HPX::wrap_main)
target_include_directories (single_cell PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/utils)
