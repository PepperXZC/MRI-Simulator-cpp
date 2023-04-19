#ifndef BLOCHSIM_H
#define BLOCHSIM_H

#include <eigen3/Eigen/Dense>
#include <cmath>

using Eigen::Matrix3d;
using Eigen::Vector3d;

struct Arg_Mats
{
    Matrix3d A;
    Vector3d B;
};

Matrix3d Rz(double theta);
Matrix3d Ry(double theta);
Matrix3d Rx(double theta);

Arg_Mats freeprecess(double T, double T1, double T2, double df);

#endif