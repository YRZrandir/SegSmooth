cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
project("SegSmooth" C CXX)

set(CMAKE_CXX_FLAGS"${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_C_FLAGS"${CMAKE_C_FLAGS} -fopenmp")

find_package(OpenMP REQUIRED)
find_package(CGAL REQUIRED)
find_package(Ceres REQUIRED)
find_package(assimp REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(SegSmooth "main.cpp" "Polyhedron.cpp" "Polyhedron.h" "SegSmooth.cpp" "SegSmooth.h")

target_compile_definitions(SegSmooth PRIVATE CGAL_PMP_USE_CERES_SOLVER)
target_link_libraries(SegSmooth CGAL::CGAL assimp::assimp nlohmann_json::nlohmann_json Ceres::ceres)
target_include_directories(SegSmooth PRIVATE "./")