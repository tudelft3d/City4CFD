cmake_minimum_required(VERSION 3.1)

project(City4CFD)

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/" ${CMAKE_MODULE_PATH})
#set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE "Debug")

#if (COMMAND cmake_policy)
#    cmake_policy(SET CMP0003 NEW)
#endif()

#==============================
# Alpha wrap support
set(ALPHA_WRAP_SUPPORT  "Include CGAL version with Alpha Wrap"  OFF)
#==============================

# BOOST
find_package(Boost 1.66 REQUIRED COMPONENTS filesystem)

# CGAL
find_package(CGAL REQUIRED QUIET COMPONENTS)
if (CGAL_FOUND)
    if (CGAL_VERSION VERSION_GREATER_EQUAL "5.1")
        message(STATUS "Found CGAL ${CGAL_VERSION}")
    elseif (CGAL_VERSION VERSION_GREATER_EQUAL "5.0")
        message(STATUS "Found CGAL version greater than 5.0, but less than 5.1. "
                "Proceeding to compile with included CGAL headers.")
        include_directories(${CMAKE_SOURCE_DIR}/thirdparty/CGAL/include)
        set(CGAL_USE_INCLUDED_HEADERS ON)
    else()
        message(FATAL_ERROR "Found CGAL version ${CGAL_VERSION} which is not supported!"
                "Please use CGAL version 5")
        return()
    endif ()
endif ()
if (NOT CGAL_USE_INCLUDED_HEADERS AND ALPHA_WRAP_SUPPORT)
    message(STATUS "Requested alpha wrap support. Compiling with included CGAL headers.")
    include_directories( ${CMAKE_SOURCE_DIR}/thirdparty/CGAL/include )
endif()

# LASlib with CGAL and LAStools
#include_directories(${CMAKE_SOURCE_DIR}/thirdparty/LAStools/LASlib/inc)
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/LAStools)
find_package(LASLIB REQUIRED)
file(GLOB LASZIP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/LAStools/LASzip/src)
file(GLOB LASLIB_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/LAStools/LASlib/inc)
if(LASLIB_FOUND AND NOT TARGET CGAL::LASLIB_support)
    add_library(CGAL::LASLIB_support INTERFACE IMPORTED)
    set_target_properties(CGAL::LASLIB_support PROPERTIES
            INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_LASLIB"
            INTERFACE_INCLUDE_DIRECTORIES "${LASLIB_INCLUDE_DIR};${LASZIP_INCLUDE_DIR}"
            INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${LASLIB_INCLUDE_DIR};${LASZIP_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "${LASLIB_LIBRARIES}")
endif()

# Other third-party
include_directories(${CMAKE_SOURCE_DIR}/thirdparty)
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/valijson)
include_directories(${Boost_INCLUDE_DIR})

# Creating entries for target: City4CFD
FILE(GLOB SRC_FILES "src/*.cpp")
add_executable(City4CFD ${SRC_FILES})
set_target_properties(
        City4CFD
        PROPERTIES CXX_STANDARD 14
)

target_link_libraries(City4CFD
        ${CGAL_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_FILESYSTEM_LIBRARY}
        LASlib
        CGAL::LASLIB_support
        )

install(TARGETS City4CFD DESTINATION bin)