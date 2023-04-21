#include "info.h"

proton_info::proton_info(bool flow, double Mx, double My,double Mz)
{
    this->flow = flow;
    this->Mxy = {Mx, My};
    this->amplitude = abs(this->Mxy);
    this->phase = arg(this->Mxy);
    this->M << Mx, My, Mz;
}

proton_info::proton_info(bool flow, Vector3d M)
{
    this->M = M;
    this->flow = flow;
    this->Mxy = {M(0), M(1)};
    this->amplitude = abs(this->Mxy);
    this->phase = arg(this->Mxy);
}

pool_info::pool_info()
{
    // 应该都是初始化为0？
}

pool_info::pool_info(double fov, double delta, double z0, double HR, double flow_speed, double bandwidth, double TR, int num_vassels)
{
    this->fov = fov;
    this->delta = delta;
    this->z0 = z0;                                            
    this->HR = HR;
    this->bandwidth = bandwidth;        // cm
    this->flow_speed = flow_speed;
    this->N_read = int(fov / delta);
    this->N_pe = int(fov / delta);

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

void pool_info::T_info_generate(double T1_vassel, double T2_vassel, double T1_tissue, double T2_tissue)
{
    T_vassel[0] = T1_vassel; T_vassel[1] = T2_vassel;
    T_tissue[0] = T1_tissue; T_tissue[1] = T2_tissue;
}

bool operator<(const Vector3d& vct1, const Vector3d& vct2)
{
    if (vct1.sum() < vct2.sum()) return true;
    else return false;
}

Voxel::Voxel()
{
    this->M << 0, 0, 1;
    this->num_protons = 1;
    this->flow = 0;
}

void Voxel::initialize(double T[], bool flow, Vector3d position, Vector3d M, int num_protons = 1, int index_vassel = 0)
{
    this->T1 = T[0];
    this->T2 = T[1];
    this->position = position;
    this->M = M;
    this->flow_speed = 0;
    this->flow = flow;
    this->num_protons = num_protons;
    proton_init(M);
}

void Voxel::initialize(const Voxel &vox)
{
    this->T1 = vox.T1;
    this->T2 = vox.T2;
    this->position = vox.position;
    this->M = vox.M;
    this->flow_speed = vox.flow_speed;
    this->flow = vox.flow;
    this->num_protons = vox.num_protons;
    proton_init(vox.M);
}

Voxel::~Voxel()
{
}

void Voxel::proton_init(Vector3d M)
{
    for (int i = 0; i < this->num_protons; i++)
    {
        proton_info proton(this->flow, M);
        this->protons.push_back(proton);
    }
}

void Voxel::proton_sum()
{
    double x = 0, y = 0, z = 0;
    for (int i = 0; i < protons.size(); i++)
    {
        x = x + protons[i].M(0);
        y = y + protons[i].M(1);
        z = z + protons[i].M(2);
    }
    this->M << x / protons.size(), y / protons.size(), z / protons.size();
}

pool::pool(const pool_info& info)
{
    pool_args = info;
    pool_length = int(info.fov / info.delta);
    // vassel : centered
    index_generate();
    data_initialize();
}

pool::~pool()
{
}

void pool::index_generate()
{
    // one vassel
    int width = pool_args.bandwidth / pool_args.delta;
    int center = int(pool_length / 2); int half = int(width / 2);     // usually even center
    int lower = center - half; int upper = center + half;
    std::cout << lower << " " << upper << std::endl;
    for (int i = 0; i < pool_length; i++){
        for (int j = 0; j < pool_length; j ++){
            for (int k = 0; k < pool_length; k++){          // z : from bottom to top
                Vector3d pos; pos << i, j, k;
                if (i >= lower && i < upper) // 27 37 longitudinal
                    this->vassel_index_vector.push_back(pos);
                else
                    this->tissue_index_vector.push_back(pos);
                // std::cout << i << " " << j << " " << k << std::endl;
            }
        }
    }
}

Vector3d pool::index_to_position(Vector3d index)
{
    Vector3d res;
    double x = - pool_args.fov / 2 + index(0) * pool_args.delta;
    double y = pool_args.fov / 2 - index(1) * pool_args.delta; // y == 0 : index == 32
    double z = index(2) * pool_args.delta;
    res << x, y, z;
    return res;
}

void pool::data_initialize()
{
    body = new Voxel **[pool_length];
    for (int i = 0; i < pool_length; i++){
        body[i] = new Voxel *[pool_length];
        for (int j = 0; j < pool_length; j++){
            body[i][j] = new Voxel [pool_length];
        }
    }
    int s = 1;
    for (int i = 0; i < vassel_index_vector.size(); i++)
    {
        Vector3d pos_index = vassel_index_vector[i];
        int x = int(pos_index(0)); int y = int(pos_index(1)); int z = int(pos_index(2));
        Vector3d M_init = {0, 0, 1}; 
        body[x][y][z].initialize(pool_args.T_vassel, 1, index_to_position(pos_index), M_init);
        // std::cout << x << " " << y << " " << z << " " << index_to_position(pos_index)(0) << " " << index_to_position(pos_index)(1) << " " << index_to_position(pos_index)(2) << " " << std::endl;
        s += 1;
    }
    for (int i = 0; i < tissue_index_vector.size(); i++)
    {
        Vector3d pos_index = tissue_index_vector[i];
        int x = int(pos_index(0)); int y = int(pos_index(1)); int z = int(pos_index(2));
        Vector3d M_init = {0, 0, 1}; 
        body[x][y][z].initialize(pool_args.T_tissue, 0, index_to_position(pos_index), M_init);
        s += 1;
    }
    std::cout << s << std::endl;
}

// friend bool operator< (const Voxel& vox1, const Voxel& vox2)
// {
//     if (vox1.index < vox2.index)
//         return true;
//     else return false; 
// }
