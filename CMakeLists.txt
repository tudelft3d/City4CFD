cmake_minimum_required(VERSION 3.9)

project(city4cfd)

#set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_BUILD_TYPE "Release")

#if (COMMAND cmake_policy)
#    cmake_policy(SET CMP0003 NEW)
#endif()

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

if (MSVC)
    add_definitions(-DNOMINMAX)
    add_definitions("/EHsc")
endif (MSVC)

# BOOST
find_package(Boost 1.66 REQUIRED COMPONENTS filesystem locale QUIET)
if (WIN32)
    FIND_PACKAGE(Boost)
    if (Boost_FOUND)
        INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
        ADD_DEFINITIONS( "-DHAS_BOOST" )
        set(Boost_USE_STATIC_LIBS   ON)
        set(Boost_ARCHITECTURE "-x64")
    endif()
endif (WIN32)

# CGAL
find_package(CGAL REQUIRED)
if (CGAL_VERSION VERSION_GREATER_EQUAL "6.0.1")
   message(STATUS "Found CGAL version ${CGAL_VERSION}")
else()
   message(FATAL_ERROR "Found CGAL version ${CGAL_VERSION} which is not supported!"
           " Please use CGAL version 6.0.1 or higher")
   return()
endif ()

# GDAL
find_package(GDAL 3.0 REQUIRED)
if (NOT GDAL_FOUND)
    message(SEND_ERROR "City4CFD requires the GDAL library, check the README for configuration")
endif()
include_directories(${GDAL_INCLUDE_DIR})

# Eigen
find_package(Eigen3 3.3.4)
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
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/roofer/src)
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/val3dity)
include_directories(${CMAKE_SOURCE_DIR}/thirdparty)
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/valijson)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/roofer/src)

# Tools
add_subdirectory(${CMAKE_SOURCE_DIR}/tools/prepare_point_cloud)

# Creating entries for target: City4CFD
FILE(GLOB SRC_FILES "src/*.cpp")
add_executable(city4cfd ${SRC_FILES})
set_target_properties(
        city4cfd
        PROPERTIES CXX_STANDARD 17
)

target_link_libraries(city4cfd
        ${CGAL_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
        ${GDAL_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_LOCALE_LIBRARY}
        LASlib
        CGAL::Eigen3_support
        roofer_api
        val3dity
        ${OpenMP_support}
        )

install(TARGETS city4cfd DESTINATION bin)
