if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0011 NEW)
endif(COMMAND cmake_policy)

cmake_minimum_required (VERSION 3.16)
project( val3dity )

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()


add_definitions(-std=c++17)

set( CMAKE_BUILD_TYPE "Release")
set( CMAKE_CXX_FLAGS "-O2" )

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

if (MSVC)
    add_definitions("/EHsc")
endif(MSVC)

# CGAL
find_package( CGAL QUIET COMPONENTS )
if (NOT CGAL_FOUND )
  message(SEND_ERROR "val3dity requires the CGAL library")
  return()  
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Eigen
find_package(Eigen3 3.1.0 QUIET)
include(CGAL_Eigen3_support)


# Boost
find_package( Boost REQUIRED COMPONENTS filesystem QUIET)

if ( NOT Boost_FOUND )
  message(SEND_ERROR "val3dity requires the Boost library")
  return()  
endif()

# GEOS
find_package( GEOS REQUIRED )
# if linux complains: 
# find_package( GEOS )
if ( NOT GEOS_FOUND )
  message(SEND_ERROR "val3dity requires the GEOS library")
endif()

include_directories( ${GEOS_INCLUDE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty )

FILE(GLOB_RECURSE THIRDPARTY thirdparty/*.cpp)
add_library(val3dity_thirdparty STATIC ${THIRDPARTY})

FILE(GLOB SRC_FILES src/*.cpp)

list(REMOVE_ITEM SRC_FILES src/main.cpp)
add_library(val3dity STATIC ${SRC_FILES})
target_include_directories(val3dity PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
set_property(TARGET val3dity_thirdparty PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET val3dity PROPERTY POSITION_INDEPENDENT_CODE ON)


#target_link_libraries(val3dity CGAL::CGAL CGAL::Eigen3_support GEOS::geos_c val3dity_thirdparty Boost::filesystem)
target_link_libraries(val3dity CGAL::CGAL CGAL::Eigen3_support ${GEOS_LIBRARY} val3dity_thirdparty Boost::filesystem)

install(TARGETS val3dity DESTINATION bin)
