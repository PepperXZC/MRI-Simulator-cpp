#ifndef ENCODING_H
#define ENCODING_H

#include <eigen3/Eigen/Dense>
#include <cmath>
#include "../Bloch/Bloch.h"
#include "../info/info.h"
#include "../SeqLoader/SeqLoader.h"
// #include "../../include/tqdm/tqdm.h"

using Eigen::Matrix;
using Eigen::Vector3d;
typedef Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> Mat;
const int MAX_length = 999;

class Simulator: public SeqLoader, pool
{
private:
    void progress(bool flow);
    void mat_initialize();
    Voxel ***sliced_body;
    int int_z0, int_tikn;
public:
    Simulator(const SeqLoader &Seq, const pool &pl);
    ~Simulator();
    void slice_select(double z0, double thickness);         // n_dim
    void RF_pulse(double fa);
    void None_operation(double t);
    void encoding(double t, double Gx, double Gy, double Gz);
    void ADC_Readout(double t, int line_index, int sample_index, double Gx, double Gy, double Gz);
    void check_readout_record(int readout_index);
    void load_seqence();
    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> test_Mat;
    Mat data;
    
    vector<Mat> result_list;
};

#endif