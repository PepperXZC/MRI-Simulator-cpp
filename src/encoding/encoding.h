#ifndef ENCODING_H
#define ENCODING_H

#include <eigen3/Eigen/Dense>
#include <cmath>
#include "../Bloch/Bloch.h"
#include "../info/info.h"
#include "../SeqLoader/SeqLoader.h"
#include "omp.h"
// #include "../../include/tqdm/tqdm.h"

using Eigen::Matrix;
using Eigen::Vector3d;
typedef Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Mat;
const int MAX_length = 999;

class Simulator 
{
private:
    void progress(bool flow);
    void mat_initialize(int n_pe, int n_read);
    // Voxel ***sliced_body;
    bool slice_selected = 0;
    int slice_lower, slice_upper; // [slice_lower, slice_upper]
public:
    Simulator(pool &pl);
    ~Simulator();
    void slice_select(pool &pl, double z0, double thickness);         // n_dim
    void RF_pulse(pool &pl, double fa);
    void None_operation(pool &pl, double t);
    void encoding(pool &pl, double t, double Gx, double Gy, double Gz);
    void ADC_Readout(pool &pl, double t, int line_index, int sample_index, double Gx, double Gy, double Gz);
    void load_seqence(const SeqLoader &Seq, pool &pl);
    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> test_Mat;
    Mat real_data, img_data;
    
    // vector<Mat> result_list;
};

#endif