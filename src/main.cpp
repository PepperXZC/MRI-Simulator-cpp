#include<iostream>
#include<eigen3/Eigen/Dense>
#include "Bloch/Bloch.h"
#include "SeqLoader/SeqLoader.h"
#include "info/info.h"
#include "encoding/encoding.h"
#include <fstream>
#include <opencv2/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/opencv.hpp>

using std::ofstream;
int main()
{
    // std::string seq_path = "./sequence/bSSFP_TR2o8.yaml";
    std::string seq_path = "/home/xzc/cpp/sequence/bSSFP_TR2o8.yaml";
    SeqLoader sequence(seq_path);


    pool_info test_pool_info(6.4, 0.1, 3.2, 60, 0, 1, 2.8, 1);
    test_pool_info.T_info_generate(1500, 50, 1000, 50);
    pool test_pool(test_pool_info);
    // std::cout << test_pool.body[0][0][0].T1 << std::endl;
    Simulator test_simulator(sequence, test_pool);
    test_simulator.load_seqence();

    std::cout << test_simulator.test_Mat << std::endl;
    ofstream fout;
    fout.open("real_mat.txt", std::ios::trunc);
    fout << test_simulator.real_data << std::endl;
    fout.close();
    fout.open("img_mat.txt", std::ios::trunc);
    fout << test_simulator.img_data << std::endl;
    fout.close();
    // cv::Mat img1;
    // cv::Mat 
    // cv::eigen2cv(test_simulator.test_Mat, img1);
    // cv::imshow("now.png", img1);
    // while(1)
    // {
    //     int key = cv::waitKey(0);
        
    //         if (key == 'q')
    //         {
    //         cv::destroyAllWindows();
    //         break;
    //         }
    // }

    // cv::Mat cvM1;
    // cv::FileStorage cvfs("data.xml", cv::FileStorage::WRITE);
    // cvfs << test_simulator.data;
    // cvfs.release();
    return 0;
}