#include "Bloch/Bloch.h"
#include "SeqLoader/SeqLoader.h"
#include "encoding/encoding.h"
#include "flow/flow.h"
#include "info/info.h"
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <vector>
// #include <opencv2/core.hpp>
// #include <opencv2/core/eigen.hpp>
// #include <opencv2/opencv.hpp>

using std::ofstream;
int main() {
  // std::string seq_path = "./sequence/bSSFP_TR2o8.yaml";
  std::string seq_path =
      "/home/xzc/MRI-Simulator-cpp/sequence/bSSFP_TR2o8.yaml";
  SeqLoader sequence(seq_path);

  pool_info test_pool_info(12.8, 0.1, 3.2, 60, 20, 1, 2.8, 1, 3);
  double T2 = 50;
  double T1_tissue = 1000;
  for (int i = 1; i < test_pool_info.num_vassels; i++) {
    double T1 = 1300 + i * 200;
    vector<double> T_vassel = {T1, T2};
    test_pool_info.get_T_vassel(T_vassel);
  }
  vector<double> T_tissue = {T1_tissue, T2};
  test_pool_info.get_T_tissue(T_tissue);
  vector<double> cli = {1.5, 5.5, 10.5};
  test_pool_info.center_generate(cli);

  pool test_pool(test_pool_info, test_pool_info.fov, test_pool_info.fov,
                 test_pool_info.fov);
  test_pool.whole_init();

  molli_sequence mlsq(sequence,
                      "/home/xzc/MRI-Simulator-cpp/sequence/molliseq.yaml");
  flow_sequence flsq(mlsq, test_pool_info);
  // flsq.print_flowseq();

  Simulator test_simulator(test_pool);
  flow_experiment fl_program(flsq.flow_molli_seq, test_pool, test_simulator);
  // flow_experiment fl_program(mlsq.molli_list, test_pool, test_simulator);
  std::string save_path = "/home/xzc/MRI-Simulator-cpp/result/flow_Mz_128/";
  fl_program.save_mat(save_path);

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