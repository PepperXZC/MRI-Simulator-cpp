PROJECT (bSSFPSeq)
aux_source_directory(. DIR_SRCS)
# include_directories($(PROJECT_SOURCE_DIR)/include)

add_subdirectory(Bloch)
add_subdirectory(info)
add_subdirectory(SeqLoader)
add_subdirectory(encoding)

add_executable(bSSFPSeq main.cpp)
target_link_libraries(bSSFPSeq Bloch)
target_link_libraries(bSSFPSeq SeqLoader)
target_link_libraries(bSSFPSeq info)
target_link_libraries(bSSFPSeq encoding)
# target_link_libraries(bSSFPSeq OpenCV)
# target_link_libraries(bSSFPSeq ${OpenCV_LIBS})
# include_directories(${OpenCV_INCLUDE_DIRS})
# target_link_libraries(bSSFPSeq ${OpenCV_LIBS})
find_package(OpenCV 4.0.0 REQUIRED)



message(STATUS "OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")

add_definitions(-DOPENCV_DISABLE_EIGEN_TENSOR_SUPPORT)