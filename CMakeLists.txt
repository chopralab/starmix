cmake_minimum_required(VERSION 3.2)

project(starmix C CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(spear REQUIRED)
find_package(lemon REQUIRED)
find_package(chemfiles REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph)

add_subdirectory(lemon_spear)
