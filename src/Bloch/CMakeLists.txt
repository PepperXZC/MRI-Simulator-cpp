project(Bloch)
aux_source_directory(. DIR_LIB_SRCS)
set(SOURCE_FILES 
    Bloch.cpp
    Bloch.h)

# add_library(BlochSim ${DIR_LIB_SRCS})
add_library(Bloch ${SOURCE_FILES})
# target_link_libraries(Bloch Eigen3::Eigen)
target_link_libraries(Bloch ${Eigen3_LIBRARY})
