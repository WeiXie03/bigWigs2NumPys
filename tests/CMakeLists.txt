add_library(
    doctest
    INTERFACE
)
target_sources(
    doctest
    INTERFACE
    ${PROJECT_SOURCE_DIR}/libs/doctest/doctest.h
)
target_include_directories(
    doctest
    INTERFACE
    ${PROJECT_SOURCE_DIR}/libs/doctest
)

add_executable(
    test
    main.cc
)

target_include_directories(test PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test libBigWig torch doctest ${PROJECT_NAME})