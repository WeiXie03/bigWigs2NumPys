file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/libbigWigs2tensors/*.h")

add_library(libbigWigs2tensors SHARED
    util.cc proc_bigWigs.cc
    ${HEADER_LIST}
)

target_include_directories(libbigWigs2tensors PUBLIC ../include ../extern/libBigWig)

target_link_libraries(libbigWigs2tensors
  ${TORCH_LIBRARIES} libBigWig
  )

target_compile_features(libbigWigs2tensors PUBLIC cxx_std_20)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})