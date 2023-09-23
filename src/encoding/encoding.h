#ifndef ENCODING_H
#define ENCODING_H

#include "../Bloch/Bloch.h"
#include "../SeqLoader/SeqLoader.h"
#include "../info/info.h"
#include "omp.h"
#include <cinttypes>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>

// #include "../../include/tqdm/tqdm.h"

using Eigen::Matrix;
using Eigen::Vector3d;
typedef Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Mat;
const int MAX_length = 999;

class Simulator {
private:
  void progress(bool flow);
  void mat_initialize(pool_info &pl, int n_pe, int n_read);
  // Voxel ***sliced_body;
  bool slice_selected = 0;

public:
  Simulator(pool &pl);
  ~Simulator();
  void slice_select(pool &sample, double z0, double thickness); // n_dim
  void RF_pulse(pool &sample, double fa);
  void None_operation(pool &sample, double t);
  void encoding(pool &sample, double t, double Gx, double Gy, double Gz);
  void ADC_Readout(pool &sample, double t, int line_index, int sample_index,
                   double Gx, double Gy, double Gz);
  void load_seqence(const SeqLoader &Seq, pool &pl);
  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> test_Mat;
  Mat real_data, img_data;
  // vector<vector<vector<double>>> Mz_data; // 目前只考虑 tk = 1
  // 的二维情形，所以在第三维直接取[slice_lower]
  Mat Mz_data;
  vector<Mat> Mz_list;
  int slice_lower, slice_upper; // [slice_lower, slice_upper]

  // vector<Mat> result_list;
};

#endif