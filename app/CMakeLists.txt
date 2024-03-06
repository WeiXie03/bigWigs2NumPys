add_executable(bigWigs2tensors main.cc)
target_compile_features(bigWigs2tensors PRIVATE cxx_std_20)

target_include_directories(bigWigs2tensors PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_directories(bigWigs2tensors PUBLIC ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(bigWigs2tensors
    libBigWig
    bigWigs2tensors_lib
    )