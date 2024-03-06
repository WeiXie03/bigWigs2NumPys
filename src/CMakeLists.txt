file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${libbigWigs2tensors_lib_SOURCE_DIR}/include/libbigWigs2tensors_lib/*.h")

add_library(bigWigs2tensors_lib STATIC
    util.cc proc_bigWigs.cc
    ${HEADER_LIST}
)

set_property(TARGET bigWigs2tensors_lib PROPERTY OUTPUT_NAME bigWigs2tensors)

target_include_directories(bigWigs2tensors_lib PUBLIC ../include ../extern/libBigWig)

target_link_libraries(bigWigs2tensors_lib
  PUBLIC ${TORCH_LIBRARIES}
  PUBLIC libBigWig
  # PRIVATE TBB::tbb
  )

target_compile_features(bigWigs2tensors_lib PUBLIC cxx_std_20)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})