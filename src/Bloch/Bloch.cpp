#include "Bloch.h"

Matrix3d Rz(double theta)
{
    Matrix3d Rz;
    Rz << cos(theta), -sin(theta), 0,
        sin(theta), cos(theta), 0,
        0, 0, 1;
    return Rz;
}

Matrix3d Ry(double theta)
{
    Matrix3d Ry;
    Ry << cos(theta), 0, sin(theta),
        0, 1, 0,
        -sin(theta), 0, cos(theta);
    return Ry;
}

Matrix3d Rx(double theta)
{
    Matrix3d Rx;
    Rx << 1, 0, 0,
        0, cos(theta), -sin(theta),
        0, sin(theta), cos(theta);
    return Rx;
}

Arg_Mats freeprecess(double T, double T1, double T2, double df)
{
    Arg_Mats args;
    double phi = 2 * M_PI * df * T / 1000;
    double E1 = exp(-T / T1);
    double E2 = exp(-T / T2);
    args.A << E2, 0, 0,
        0, E2, 0,
        0, 0, E1;
    args.A = args.A * Rz(phi);
    args.B << 0, 0, (1 - E1);
    return args;
}