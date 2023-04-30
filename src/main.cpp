#include "Bloch/Bloch.h"
#include "SeqLoader/SeqLoader.h"
#include "encoding/encoding.h"
#include "flow/flow.h"
#include "info/info.h"
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/opencv.hpp>

using std::ofstream;
int main() {
  // std::string seq_path = "./sequence/bSSFP_TR2o8.yaml";
  std::string seq_path =
      "/home/xzc/MRI-Simulator-cpp/sequence/bSSFP_TR2o8.yaml";
  SeqLoader sequence(seq_path);

  pool_info test_pool_info(6.4, 0.1, 3.2, 60, 20, 1, 2.8, 1);
  test_pool_info.T_info_generate(1500, 50, 1000, 50);
  pool test_pool(test_pool_info, test_pool_info.fov, test_pool_info.fov,
                 test_pool_info.fov);
  test_pool.whole_init();

  molli_sequence mlsq(sequence,
                      "/home/xzc/MRI-Simulator-cpp/sequence/molliseq.yaml");
  flow_sequence flsq(mlsq, test_pool_info);
  // flsq.print_flowseq();

  Simulator test_simulator(test_pool);
  flow_experiment fl_program(flsq, test_pool, test_simulator);

  // test_simulator.load_seqence(sequence, test_pool);

  // std::cout << test_simulator.test_Mat << std::endl;
  // ofstream fout;
  // fout.open("real_mat.txt", std::ios::trunc);
  // fout << test_simulator.real_data << std::endl;
  // fout.close();
  // fout.open("img_mat.txt", std::ios::trunc);
  // fout << test_simulator.img_data << std::endl;
  // fout.close();
  return 0;
}