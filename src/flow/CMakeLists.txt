project(flow)
aux_source_directory(. DIR_LIB_SRCS)
set(SOURCE_FILES
    flow.h
    flow.cpp)

add_library(flow ${SOURCE_FILES})
# target_link_libraries(Bloch Eigen3::Eigen)
target_link_libraries(flow ${Eigen3_LIBRARY})
target_link_libraries(flow info)
target_link_libraries(flow encoding)
target_link_libraries(flow SeqLoader)
target_link_libraries(flow yaml-cpp)