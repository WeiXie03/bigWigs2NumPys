cmake_minimum_required(VERSION 3.15)
project(bigWigs2tensors VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)

# include(CTest)
# enable_testing()

find_package(Torch REQUIRED)
message(STATUS "Torch found: ${TORCH_FOUND}")
message(STATUS "Torch_INCLUDE_DIRS: ${TORCH_INCLUDE_DIRS}")

include_directories(${TORCH_INCLUDE_DIRS})
add_subdirectory(libs/libBigWig)

add_library(${PROJECT_NAME} SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/src/proc_bigWigs.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/bigWigs2tensors/proc_bigWigs.h
)

# add_executable(bigWigs2tensors main.cc)

add_subdirectory(tests)
# enable_testing()
add_test(NAME test COMMAND test)

# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)
