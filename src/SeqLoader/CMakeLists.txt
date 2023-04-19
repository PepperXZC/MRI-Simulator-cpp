project(SeqLoader)

aux_source_directory(. DIR_LIB_SRCS)
set(SOURCE_FILES 
    SeqLoader.cpp
    SeqLoader.h)

add_library(SeqLoader ${SOURCE_FILES})
# target_link_libraries(Bloch Eigen3::Eigen)
target_link_libraries(SeqLoader ${Eigen3_LIBRARY})
target_link_libraries(SeqLoader yaml-cpp)