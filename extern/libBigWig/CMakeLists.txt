add_library(libBigWig SHARED IMPORTED GLOBAL)

set_target_properties(libBigWig PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/lib/libBigWig.so"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(libBigWig INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")