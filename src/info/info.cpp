#include "info.h"
#include <eigen3/Eigen/src/Core/Matrix.h>

proton_info::proton_info(bool flow, double Mx, double My, double Mz) {
  this->flow = flow;
  this->Mxy = {Mx, My};
  this->amplitude = abs(this->Mxy);
  this->phase = arg(this->Mxy);
  this->M << Mx, My, Mz;
}

proton_info::proton_info(bool flow, Vector3d M) {
  this->M = M;
  this->flow = flow;
  this->Mxy = {M(0), M(1)};
  this->amplitude = abs(this->Mxy);
  this->phase = arg(this->Mxy);
}

pool_info::pool_info() {
  // 应该都是初始化为0？
}

pool_info::pool_info(double fov, double delta, double z0, double HR,
                     double flow_speed, double bandwidth, double TR,
                     int num_protons, int num_vassels) {
  this->fov = fov;
  this->delta = delta;
  this->z0 = z0;
  this->HR = HR;
  this->bandwidth = bandwidth; // cm
  this->flow_speed = flow_speed;
  this->N_read = int(fov / delta);
  this->N_pe = int(fov / delta);
  this->num_vassels = num_vassels;
  this->num_protons = num_protons;

  // bSSFP
  this->tau_y = TR / 4;
  this->tau_x = TR / 2;
  this->delta_t = this->tau_x / this->N_read;
  this->G_x = 1e3 / (g * this->delta_t * fov);
  this->delta_ky = 1 / fov;
  this->ky_max = 0.5 * this->N_pe * this->delta_ky;
  this->G_yp = 1e3 * this->ky_max / (g * this->tau_y);
  this->G_yi = 1e3 * this->delta_ky / (g * this->tau_y);
}

void pool_info::get_T_vassel(const vector<double> &T) { T_vassel.push_back(T); }

void pool_info::get_T_tissue(const vector<double> &T) { T_tissue = T; }

void pool_info::center_generate(const vector<double> &center_li) {
  for (int i = 0; i < center_li.size(); i++)
    center_list.push_back(center_li[i]);
}

bool operator<(const Vector3d &vct1, const Vector3d &vct2) {
  if (vct1.sum() < vct2.sum())
    return true;
  else
    return false;
}

Voxel::Voxel() {
  this->M << 0, 0, 1;
  this->num_protons = 1;
  this->flow = 0;
}

void Voxel::initialize(double T1, double T2, bool flow, Vector3d position,
                       Vector3d M, int num_protons = 1) {
  this->T1 = T1;
  this->T2 = T2;
  this->position = position;
  this->M = M;
  this->flow_speed = 0;
  this->flow = flow;
  this->num_protons = num_protons;
  proton_init(M);
  proton_sum();
}

void Voxel::initialize(const Voxel &vox) {
  this->T1 = vox.T1;
  this->T2 = vox.T2;
  this->position = vox.position;
  this->M = vox.M;
  this->flow_speed = vox.flow_speed;
  this->flow = vox.flow;
  this->num_protons = vox.num_protons;
  proton_init(vox.M);
  proton_sum();
}

Voxel::~Voxel() {}

void Voxel::proton_init(Vector3d M) {
  this->protons.clear();
  for (int i = 0; i < this->num_protons; i++) {
    proton_info proton(this->flow, M);
    this->protons.push_back(proton);
  }
}

void Voxel::proton_sum() {
  double x = 0, y = 0, z = 0;
  for (int i = 0; i < protons.size(); i++) {
    x = x + protons[i].M(0);
    y = y + protons[i].M(1);
    z = z + protons[i].M(2);
  }
  this->M << x / protons.size(), y / protons.size(), z / protons.size();
}

pool::pool(const pool_info &info, double x, double y, double z) {
  pool_args = info;
  // pool_length = int(info.fov / info.delta);
  x_length = x / info.delta;
  y_length = y / info.delta;
  z_length = z / info.delta;
  // empty_init();
  // vassel : centered
}

void pool::vassel_init(double T1, double T2, int lower, double y_pos) {
  body = new Voxel **[x_length];
  for (int i = 0; i < x_length; i++) {
    body[i] = new Voxel *[y_length];
    for (int j = 0; j < y_length; j++) {
      body[i][j] = new Voxel[z_length];
      for (int k = 0; k < z_length; k++) {
        Vector3d M = {0, 0, 1};
        Vector3d pos = {-pool_args.fov / 2 + (i + lower) * pool_args.delta,
                        y_pos, k * pool_args.delta};
        body[i][j][k].initialize(T1, T2, 1, pos, M, pool_args.num_protons);
      }
    }
  }
}

pool::~pool() {}

void pool::get_lower_upper(const vector<double> &center_list) {
  for (double c : center_list) {
    int _center = c / pool_args.delta;
    int _half = pool_args.bandwidth / (2 * pool_args.delta);
    // int half = int(width / 2); // usually even center
    lower_list.push_back(_center - _half);
    upper_list.push_back(_center + _half);
    std::cout << _center - _half << " " << _center + _half << std::endl;
  }
}

void pool::index_generate(const vector<double> &center_list) {
  // one vassel
  // int width = pool_args.bandwidth / pool_args.delta;
  // int center = x_length / 2;
  vassel_list_init();
  get_lower_upper(center_list);
#pragma omp parallel for
  for (int i = 0; i < x_length; i++) {
    for (int j = 0; j < y_length; j++) {
      for (int k = 0; k < z_length; k++) { // z : from bottom to top
        Vector3d _pos;
        _pos << i, j, k;
        bool flag = 0;
        for (int index = 0; index < pool_args.num_vassels; index++) {
          if (i >= lower_list[index] && i < upper_list[index]) {
            vassel_index_list[index].push_back(_pos);
            flag = 1;
            break;
          } // [59 69) longitudinal
        }
        if (flag == 0)
          tissue_index_vector.push_back(_pos);
      }
    }
  }
}

void pool::vassel_list_init() {
  vassel_index_list.resize(pool_args.num_vassels);
}

bool pool::check_vassel(pool &pl, double i, double j, double k) {
  int flag = 0;
  // in-plane pattern
  for (int k = 0; k < vassel_index_list.size(); k++) {
    if (i >= lower_list[k] && i < upper_list[k])
      flag = 1;
  }
  if (flag) {
    return true;
  } else
    return false;
}

void pool::whole_init() {
  index_generate(pool_args.center_list);
  data_initialize();
}

Vector3d pool::index_to_position(Vector3d index) {
  Vector3d res;
  double x = -pool_args.fov / 2 + index(0) * pool_args.delta;
  double y =
      pool_args.fov / 2 - index(1) * pool_args.delta; // y == 0 : index == 32
  double z = index(2) * pool_args.delta;
  res << x, y, z;
  return res;
}

void pool::data_initialize() {
  body = new Voxel **[x_length];
  for (int i = 0; i < x_length; i++) {
    body[i] = new Voxel *[y_length];
    for (int j = 0; j < y_length; j++) {
      body[i][j] = new Voxel[z_length];
    }
  }
  for (int k = 0; k < pool_args.center_list.size(); k++) {
    for (int i = 0; i < vassel_index_list[k].size(); i++) {
      Vector3d pos_index = vassel_index_list[k][i];
      int x = pos_index(0);
      int y = pos_index(1);
      int z = pos_index(2);
      Vector3d M_init = {0, 0, 1};
      body[x][y][z].initialize(
          pool_args.T_vassel[k][0], pool_args.T_vassel[k][1], 1,
          index_to_position(pos_index), M_init, pool_args.num_protons);
      // std::cout << x << " " << y << " " << z << " " <<
      // index_to_position(pos_index)(0) << " " <<
      // index_to_position(pos_index)(1)
      // << " " << index_to_position(pos_index)(2) << " " << std::endl;
      // s += 1;
    }
  }
  for (int i = 0; i < tissue_index_vector.size(); i++) {
    Vector3d pos_index = tissue_index_vector[i];
    int x = pos_index(0);
    int y = pos_index(1);
    int z = pos_index(2);
    Vector3d M_init = {0, 0, 1};
    body[x][y][z].initialize(pool_args.T_tissue[0], pool_args.T_tissue[1], 0,
                             index_to_position(pos_index), M_init,
                             pool_args.num_protons);
  }
}

void pool::pool_roll(int index_vassel) {
  // std::cout << "before flow: " << body[30][61][32].M(1) << " "
  //           << body[30][62][32].M(1) << " " << body[30][63][32].M(1);
  // std::cout << "pos: " << body[30][61][32].position << " "
  //           << body[30][62][32].position << " " << body[30][63][32].position
  //           << std::endl;
  for (int i = lower_list[index_vassel]; i < upper_list[index_vassel]; i++)
    for (int j = y_length - 1; j > 0; j--)
      for (int k = 0; k < z_length; k++)
        body[i][j][k].initialize(pool_args.T_vassel[index_vassel][0],
                                 pool_args.T_vassel[index_vassel][1], 1,
                                 body[i][j][k].position, body[i][j - 1][k].M,
                                 1);
  // std::cout << "after flow: " << body[30][61][32].M(1) << " "
  //           << body[30][62][32].M(1) << " " << body[30][63][32].M(1);
  // std::cout << "pos: " << body[30][61][32].position << " "
  //           << body[30][62][32].position << " " << body[30][63][32].position
  //           << std::endl;
}
