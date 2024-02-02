add_executable(bigWigs2tensors main.cc)

target_include_directories(bigWigs2tensors PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_directories(bigWigs2tensors PUBLIC ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(bigWigs2tensors
    libBigWig
    libbigWigs2tensors
    )