cmake_minimum_required(VERSION 3.10)
project(CSF LANGUAGES CXX)

find_package(OpenMP)

add_subdirectory(src)

set(BUILD_SHARED_LIBS "Build as shared library" OFF)