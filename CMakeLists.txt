cmake_minimum_required(VERSION 3.16)

project(genome_mapper LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(mapper main.cpp)

# Optimizations
target_compile_options(mapper PRIVATE -O3)

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(mapper PRIVATE OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP not found. Building without OpenMP; performance may be lower.")
endif()
