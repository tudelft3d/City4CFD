project(CityCFD)

cmake_minimum_required(VERSION 3.1)

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/" ${CMAKE_MODULE_PATH})
#set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

set( CMAKE_CXX_FLAGS "-O3" )
set( CMAKE_BUILD_TYPE "Release" )
#set( CMAKE_BUILD_TYPE "Debug" )

#if ( COMMAND cmake_policy )
#    cmake_policy( SET CMP0003 NEW )
#endif()

# BOOST
find_package( Boost REQUIRED COMPONENTS filesystem)

# LASLIB
#find_package( LASlib REQUIRED )

# CGAL
find_package( CGAL REQUIRED QUIET COMPONENTS) #remove qt5

# Third-party
include_directories( ${CMAKE_SOURCE_DIR}/thirdparty )
include_directories( ${CMAKE_SOURCE_DIR}/thirdparty/valijson )

# Creating entries for target: cityCFD
FILE(GLOB SRC_FILES "src/*.cpp")
add_executable(CityCFD ${SRC_FILES})
set_target_properties(
        CityCFD
        PROPERTIES CXX_STANDARD 14
)

target_link_libraries( CityCFD ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} Boost::filesystem)

install(TARGETS CityCFD DESTINATION bin)