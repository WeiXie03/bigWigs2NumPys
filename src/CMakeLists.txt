file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${libbigWigs2tensors_SOURCE_DIR}/include/libbigWigs2tensors/*.h")

add_library(libbigWigs2tensors STATIC
    util.cc proc_bigWigs.cc
    ${HEADER_LIST}
)

target_include_directories(libbigWigs2tensors PUBLIC ../include ../extern/libBigWig)

target_link_libraries(libbigWigs2tensors
  PUBLIC ${TORCH_LIBRARIES}
  PUBLIC libBigWig
  PRIVATE TBB::tbb
  )

target_compile_features(libbigWigs2tensors PUBLIC cxx_std_20)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})