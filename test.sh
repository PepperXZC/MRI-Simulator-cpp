#!/bin/bash
cd ./build
cmake ..
make 
# g++ ../src/main.cpp -o cv -I/usr/include/opencv4
sh ~/cpp/build/bin/bSSFPSeq
