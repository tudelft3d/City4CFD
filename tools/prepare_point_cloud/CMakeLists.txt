cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_BUILD_TYPE "Release")

#if (COMMAND cmake_policy)
#    cmake_policy(SET CMP0003 NEW)
#endif()

# OpenMP
find_package(OpenMP)

# Creating entries for target: City4CFD
add_executable(city4cfd_prepare_point_cloud src/city4cfd_prepare_point_cloud.cpp)

target_link_libraries(city4cfd_prepare_point_cloud
        LASlib
        CSF
        OpenMP::OpenMP_CXX
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_LOCALE_LIBRARY}
        )

install(TARGETS city4cfd_prepare_point_cloud DESTINATION bin)
