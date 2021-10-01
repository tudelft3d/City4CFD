project(cityCFD)

cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)

#set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

set( CMAKE_CXX_FLAGS "-O3" )
set( CMAKE_BUILD_TYPE "Release" )
#set( CMAKE_BUILD_TYPE "Debug" )

#if ( COMMAND cmake_policy )
#    cmake_policy( SET CMP0003 NEW )
#endif()

# LASLIB
find_package( LASlib REQUIRED )

# CGAL
find_package( CGAL REQUIRED QUIET COMPONENTS )

include_directories( ${CMAKE_SOURCE_DIR}/thirdparty )

# Creating entries for target: cityCFD
FILE(GLOB SRC_FILES "src/*.cpp")
add_executable(cityCFD ${SRC_FILES})
set_target_properties(
        cityCFD
        PROPERTIES CXX_STANDARD 17
)

#target_link_libraries( cityCFD ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${GDAL_LIBRARY} ${YAML_CPP_LIBRARIES} Boost::program_options Boost::filesystem Boost::locale Boost::chrono LASlib)
target_link_libraries( cityCFD ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

install(TARGETS cityCFD DESTINATION bin)