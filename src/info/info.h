#ifndef INFO_H
#define INFO_H
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <map>
#include <vector>
#include <complex>
#include <cmath>
#define g 4258 // Hz/G

using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::VectorXd;
// using Eigen::Matrix;

using std::complex;
using std::vector;
using std::map;
const int N = 99;
    
typedef struct T_info
{
    double T1;
    double T2;
}TI_info;

// 每个质子的info
struct proton_info
{
    bool flow;
    Vector3d M;
    double Mz;
    complex<double> Mxy;
    double phase;
    double amplitude;
    proton_info(bool flow, double Mx, double My,double Mz);
    proton_info(bool flow, Vector3d M);
};

struct pool_info
{
    double fov;
    double delta;           // spatial resolution
    double z0;              // slice center: cm      
    double delta_k;
    double HR;
    double B0;
    double bandwidth;       // vassel width : cm
    double flow_speed;
    int N_read;
    int N_pe;
    int num_vassels;
    // Matrix<double, Eigen::Dynamic, 2> T_vassel;
    double T_vassel[2];
    double T_tissue[2];
    // vector<double>T_tissue(2);

    double tau_y;
    double tau_x;
    double delta_t;
    double G_x;
    double delta_ky;
    double ky_max;
    double G_yp;
    double G_yi;
    // pool_info(double fov, double delta, double z0,double thickness, double HR, double bandwidth, double TR, int N_read, int N_pe); // bSSFP: determined by TR
    pool_info();
    pool_info(double fov, double delta, double z0, double HR, double flow_speed, double bandwidth, double TR, int num_vassels);
    void T_info_generate(double T1_vassel, double T2_vassel, double T1_tissue, double T2_tissue);
};

bool operator<(const Vector3d& vct1, const Vector3d& vct2);

class Voxel
{
private:
    void proton_init(Vector3d M);
public:
    double T1;
    double T2;
    Vector3d position;      // real position in fov
    Vector3d M;
    bool flow;
    double flow_speed;
    int num_protons;            // only for <
    vector<proton_info> protons;
    Voxel();
    void initialize(double T[], bool flow, Vector3d position, Vector3d M, int num_protons, int index_vassel);
    void initialize(const Voxel &vox);
    void proton_sum();
    ~Voxel();
};

class pool
{
private:
    
    // Matrix<vector<Voxel>, Eigen::Dynamic, Eigen::Dynamic> body;
    // Voxel*** body;
    void index_generate();
    void data_initialize();
public:
    int pool_length;
    vector<Vector3d> vassel_index_vector, tissue_index_vector;
    Vector3d index_to_position(Vector3d index);
    
    pool_info pool_args;
    Voxel ***body;
    pool(const pool_info& info);
    ~pool();
};
#endif