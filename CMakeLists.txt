cmake_minimum_required(VERSION 3.1)

project(city4cfd)

#set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_BUILD_TYPE "Release")

#if (COMMAND cmake_policy)
#    cmake_policy(SET CMP0003 NEW)
#endif()

# BOOST
find_package(Boost 1.66 REQUIRED COMPONENTS filesystem locale)

# CGAL
find_package(CGAL REQUIRED QUIET COMPONENTS)
if (CGAL_VERSION VERSION_GREATER_EQUAL "5.0")
    message(STATUS "Found CGAL version greater than 5.0")
    include_directories(${CMAKE_SOURCE_DIR}/thirdparty/CGAL/include)
    set(CGAL_USE_INCLUDED_HEADERS ON)
else()
    message(FATAL_ERROR "Found CGAL version ${CGAL_VERSION} which is not supported!"
            "Please use CGAL version 5")
    return()
endif ()

# Eigen
find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
include(${CMAKE_SOURCE_DIR}/thirdparty/CGAL/cmake/CGAL_Eigen3_support.cmake)

# OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    set(OpenMP_support "OpenMP::OpenMP_CXX")
else()
    message(STATUS "OpenMP not found. All functions will run single-threaded")
endif()

# Third-party
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/LAStools)
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/CSF/src)
include_directories(${CMAKE_SOURCE_DIR}/thirdparty)
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/valijson)
include_directories(${Boost_INCLUDE_DIR})

# Tools
add_subdirectory(${CMAKE_SOURCE_DIR}/tools/prepare_point_cloud)

# Creating entries for target: City4CFD
FILE(GLOB SRC_FILES "src/*.cpp")
add_executable(city4cfd ${SRC_FILES})
set_target_properties(
        city4cfd
        PROPERTIES CXX_STANDARD 14
)

target_link_libraries(city4cfd
        ${CGAL_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_LOCALE_LIBRARY}
        LASlib
        CGAL::Eigen3_support
        ${OpenMP_support}
        )

install(TARGETS city4cfd DESTINATION bin)