#ifndef INFO_H
#define INFO_H
#include <cmath>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <map>
#include <vector>
#define g 4258 // Hz/G

using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::VectorXd;
// using Eigen::Matrix;

using std::complex;
using std::map;
using std::vector;
const int N = 99;

typedef struct T_info {
  double T1;
  double T2;
} TI_info;

// 每个质子的info
struct proton_info {
  bool flow;
  Vector3d M;
  double Mz;
  complex<double> Mxy;
  double phase;
  double amplitude;
  proton_info(bool flow, double Mx, double My, double Mz);
  proton_info(bool flow, Vector3d M);
};

struct pool_info {
  double fov;
  double delta; // spatial resolution
  double z0;    // slice center: cm
  double delta_k;
  double HR;
  double B0;
  double bandwidth; // vassel width : cm
  double flow_speed;
  int N_read;
  int N_pe;
  int num_vassels;
  int num_protons;
  // Matrix<double, Eigen::Dynamic, 2> T_vassel;
  vector<vector<double>> T_vassel; // [T1, T2, ...]
  vector<double> T_tissue;         // [T1, T2, ...]
  // vector<double>T_tissue(2);

  double tau_y;
  double tau_x;
  double delta_t;
  double G_x;
  double delta_ky;
  double ky_max;
  double G_yp;
  double G_yi;
  vector<double> center_list;
  // pool_info(double fov, double delta, double z0,double thickness, double HR,
  // double bandwidth, double TR, int N_read, int N_pe); // bSSFP: determined by
  // TR
  pool_info();
  pool_info(double fov, double delta, double z0, double HR, double flow_speed,
            double bandwidth, double TR, int num_protons, int num_vassels);
  void get_T_vassel(const vector<double> &T);
  void get_T_tissue(const vector<double> &T);
  void center_generate(const vector<double> &center_li);
};

bool operator<(const Vector3d &vct1, const Vector3d &vct2);

class Voxel {
private:
  void proton_init(Vector3d M);

public:
  double T1;
  double T2;
  Vector3d position; // real position in fov
  Vector3d M;
  bool flow;
  double flow_speed;
  int num_protons; // only for <
  vector<proton_info> protons;
  Voxel();
  void initialize(double T1, double T2, bool flow, Vector3d position,
                  Vector3d M, int num_protons);
  void initialize(const Voxel &vox);
  void proton_sum();
  ~Voxel();
};

class pool {
private:
  // Matrix<vector<Voxel>, Eigen::Dynamic, Eigen::Dynamic> body;
  // Voxel*** body;
  void index_generate(const vector<double> &center); // center, bandwidth: cm
  void data_initialize();

public:
  vector<int> lower_list, upper_list;
  int x_length, y_length, z_length;
  // int lower, upper; // only one vassel
  //   int pool_length;
  // vector<Vector3d> vassel_index_vector, tissue_index_vector;
  vector<vector<Vector3d>> vassel_index_list;
  vector<Vector3d> tissue_index_vector;
  Vector3d index_to_position(Vector3d index);
  bool check_vassel(pool &pl, double i, double j, double k);
  void get_lower_upper(const vector<double> &center_list);
  void vassel_list_init();
  void vassel_init(double T1, double T2, int lower, double y_pos);
  void whole_init();
  void pool_roll(int index_vassel);
  pool_info pool_args;
  Voxel ***body = NULL;
  pool(const pool_info &info, double x, double y, double z);
  ~pool();
};
#endif