cmake_minimum_required(VERSION 3.2)

project(starmix C CXX)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

find_package(spear REQUIRED)
find_package(lemon REQUIRED)
find_package(chemfiles REQUIRED)
find_package(OpenMM REQUIRED)
find_package(pugixml REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph)

add_subdirectory(lemon_spear)
