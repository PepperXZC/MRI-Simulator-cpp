project(info)
aux_source_directory(. DIR_LIB_SRCS)
set(SOURCE_FILES 
    info.cpp
    info.h)

# add_library(BlochSim ${DIR_LIB_SRCS})
add_library(info ${SOURCE_FILES})
# target_link_libraries(Bloch Eigen3::Eigen)
target_link_libraries(info ${Eigen3_LIBRARY})