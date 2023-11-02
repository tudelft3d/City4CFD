cmake_minimum_required(VERSION 3.9)
project(city4cfd_pcprep)

set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_BUILD_TYPE "Release")

#if (COMMAND cmake_policy)
#    cmake_policy(SET CMP0003 NEW)
#endif()

# Creating entries for target: City4CFD
add_executable(city4cfd_pcprep src/city4cfd_pcprep.cpp)
set_target_properties(
        city4cfd_pcprep
        PROPERTIES CXX_STANDARD 14
)

target_link_libraries(city4cfd_pcprep
        LASlib
        CSF
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_LOCALE_LIBRARY}
        )

install(TARGETS city4cfd_pcprep DESTINATION bin)
