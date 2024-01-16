file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/lib${PROJECT_NAME}/*.h")

add_library(lib${PROJECT_NAME} SHARED
    util.cc proc_bigWigs.cc
    ${HEADER_LIST}
)

target_include_directories(lib${PROJECT_NAME} PUBLIC ../include ../extern/libBigWig)

target_link_libraries(lib${PROJECT_NAME}
  ${TORCH_LIBRARIES} libBigWig
  )

target_compile_features(lib${PROJECT_NAME} PUBLIC cxx_std_20)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})