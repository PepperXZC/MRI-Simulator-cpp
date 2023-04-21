project(encoding)
aux_source_directory(. DIR_LIB_SRCS)
set(SOURCE_FILES 
    encoding.cpp
    encoding.h)

# add_library(BlochSim ${DIR_LIB_SRCS})
add_library(encoding ${SOURCE_FILES})

include_directories($(PROJECT_SOURCE_DIR)/include)
# target_link_libraries(Bloch Eigen3::Eigen)
target_link_libraries(encoding ${Eigen3_LIBRARY})
target_link_libraries(encoding Bloch)
target_link_libraries(encoding info)
target_link_libraries(encoding SeqLoader)
