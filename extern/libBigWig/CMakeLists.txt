include(FindCURL)
include(FindZLIB)
find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)
include_directories(${CURL_INCLUDE_DIR})
include_directories(${ZLIB_INCLUDE_DIR})

# add_library(libBigWig SHARED IMPORTED GLOBAL)

# set_target_properties(libBigWig PROPERTIES
#     IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/lib/libBigWig.so"
#     INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include"
# )
add_library(libBigWig STATIC IMPORTED GLOBAL)

target_link_libraries(libBigWig INTERFACE ${CURL_LIBRARIES} ${ZLIB_LIBRARIES})

set_target_properties(libBigWig PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/lib/libBigWig.a"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(libBigWig INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")