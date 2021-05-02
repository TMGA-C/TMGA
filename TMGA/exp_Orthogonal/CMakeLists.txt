cmake_minimum_required(VERSION 3.9)
project(WFGA)

include_directories(${PROJECT_SOURCE_DIR}/include)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "-fopenmp")
add_subdirectory(${PROJECT_SOURCE_DIR}/lib)
add_executable(WFGA main.cpp include/config.h include/GenerateAChrom.h TMGA.cpp TMGA.h)
target_link_libraries(WFGA tools GenOperator pugixml classAndVarDefine config)