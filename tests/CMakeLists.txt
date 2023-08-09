add_library(
    doctest
    INTERFACE
)
target_sources(
    doctest
    INTERFACE
    ${PROJECT_SOURCE_DIR}/include/doctest/doctest.h
)
target_include_directories(
    doctest
    INTERFACE
    ${PROJECT_SOURCE_DIR}/include/doctest
)

add_executable(
    test
    main.cc
    ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}/util.h
    ${PROJECT_SOURCE_DIR}/src/util.cc
)

target_include_directories(test PUBLIC ${CMAKE_SOURCE_DIR}/include ../extern/libBigWig)
target_link_libraries(test libBigWig torch doctest lib${PROJECT_NAME})