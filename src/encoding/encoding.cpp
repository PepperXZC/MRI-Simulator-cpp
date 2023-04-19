#include "encoding.h"

void Simulator::mat_initialize()
{
    data.resize(pool_args.N_pe, pool_args.N_read);
    test_Mat.resize(pool_args.N_pe, pool_args.N_read);
}

Simulator::Simulator(const SeqLoader &Seq, const pool &pl)
    :SeqLoader(Seq), pool(pl)
{
    mat_initialize();
    // std::cout << " " <<  body[0][0][0].T1 << " " << body[0][0][0].T2 << std::endl;
}
Simulator::~Simulator()
{
}

void Simulator::slice_select(double z0, double thickness)
{
    sliced_body = new Voxel **[pool_length];
    int_z0 = z0 / pool_args.delta; int_tikn = thickness / pool_args.delta;
    for (int i = 0; i < pool_length; i++){
        sliced_body[i] = new Voxel *[pool_length];
        for (int j = 0; j < pool_length; j++){
            sliced_body[i][j] = new Voxel [int_tikn];
        }
    }
    // 30 34 32 5
    int lower = int_z0 - int_tikn / 2; int upper = int_z0 + int_tikn / 2; // n_dim
    std::cout << lower << " " << upper << " " << int_z0 <<  " " << int_tikn << std::endl;
    for (int i = 0; i < pool_length; i++){
        for (int j = 0; j < pool_length; j++){
            for (int u = 0, k = lower; k <= upper; u++, k++){
                sliced_body[i][j][u].initialize(body[i][j][k]);
                // std::cout << " " <<  body[0][0][0].T1 << " " << body[0][0][0].T2 << std::endl;
                // std::cout << " " <<  body[i][j][k].T1 << " " << body[i][j][k].T2 << std::endl;
                
            }
        }
    }
}

void Simulator::RF_pulse(double fa)
{
    // default x_rot
    Matrix3d fa_mat = Rx(fa * M_PI / 180);
    for (int i = 0; i < pool_length; i++)
        for (int j = 0; j < pool_length; j++)
            for (int k = 0; k < pool_length; k++){
                for (int u = 0; u < body[i][j][k].protons.size(); u++){
                    proton_info new_proton(body[i][j][k].flow, fa_mat * body[i][j][k].protons[u].M);
                    body[i][j][k].protons[u] = new_proton;
                } 
                body[i][j][k].proton_sum();
            }                
}

void Simulator::None_operation(double t)
{
    encoding(t, 0, 0, 0);
}

void Simulator::encoding(double t, double Gx = 0, double Gy = 0, double Gz = 0)
{
    for (int i = 0; i < pool_length; i++){
        for (int j = 0; j < pool_length; j++){
            for (int k = 0; k < int_tikn; k++){
                Vector3d pos = sliced_body[i][j][k].position;
                double gradient = g * (Gx * pos(0) + Gy * pos(1) + Gz * pos(2));
                Arg_Mats args = freeprecess(t, sliced_body[i][j][k].T1, sliced_body[i][j][k].T2, gradient);
                for (int u = 0; u < sliced_body[i][j][k].protons.size(); u++){
                    proton_info new_proton(sliced_body[i][j][k].flow, args.A * sliced_body[i][j][k].protons[u].M + args.B);
                    sliced_body[i][j][k].protons[u] = new_proton;
                }
                sliced_body[i][j][k].proton_sum();
            }
        }
    }
}

void Simulator::ADC_Readout(double t, int line_index, int sample_index, double Gx = 0, double Gy = 0, double Gz = 0)
{
    // usually readout after delta_t

    // default: center - started
    double sum_test_data = 0;
    complex<double> sum_data = {0, 0};
    for (int i = 0; i < pool_length; i++){
        for (int j = 0; j < pool_length; j++){
            for (int k = 0; k < int_tikn; k++){
                Vector3d pos = sliced_body[i][j][k].position;
                double gradient = g * (Gx * pos(0) + Gy * pos(1) + Gz * pos(2));
                Arg_Mats args = freeprecess(t, sliced_body[i][j][k].T1, sliced_body[i][j][k].T2, gradient);
                // std::cout << pos(0)  << " " << pos(1)  << " " << pos(2)  << " " << gradient << " " <<  sliced_body[i][j][k].T1 << std::endl;
                // std::cout << Gx  << " " << Gy  << " " << Gz  << std::endl;
                
                for (int u = 0; u < sliced_body[i][j][k].protons.size(); u++){
                    proton_info new_proton(sliced_body[i][j][k].flow, args.A * sliced_body[i][j][k].protons[u].M + args.B);
                    sliced_body[i][j][k].protons[u] = new_proton;
                    // std::cout << args.A * sliced_body[i][j][k].protons[u].M + args.B << std::endl;
                }
                    
                sliced_body[i][j][k].proton_sum();
                proton_info sum(sliced_body[i][j][k].flow, sliced_body[i][j][k].M);
                complex<double> Mxy = sum.Mxy;
                sum_data = sum_data + Mxy;       // only sum the mean of all protons in each voxel
                sum_test_data = sum_test_data + sum.amplitude;
            }
        }
    }
            
    data(line_index, sample_index) = sum_data;
    test_Mat(line_index, sample_index) = sum_test_data;
}

void Simulator::check_readout_record(int readout_index)
{
    if (readout_index == pool_args.N_read - 1)
    {
        result_list.push_back(data);
        mat_initialize();
    }
}

void Simulator::load_seqence()
{
    // for (int i = 0; i < tqdm::range(SeqList.size()); i++)
    for (int i = 0; i < SeqList.size(); i++)
    {
        switch (SeqList[i].type)
        {
            case NONE: {
                None_operation(SeqList[i].t);
                break;
            }
            case PULSE:{
                RF_pulse(SeqList[i].FA);
                // std::cout << "step:" << i + 1 << " " << SeqList[i].type << std::endl;
                std::cout << "PULSE!" << std::endl;
                if (i == 0)
                    slice_select(pool_args.z0, SeqList[i].thickness);
                // std::cout << " " <<  sliced_body[0][0][0].M(0) << " " << sliced_body[1][0][0].M(0) << " " << sliced_body[2][0][0].M(0) << std::endl;
                break;
            }
            case ENCODING:{
                encoding(SeqList[i].t, SeqList[i].Gx, SeqList[i].Gy, 0);
                // std::cout << " " <<  sliced_body[0][0][0].M(0) << " " << sliced_body[1][0][0].M(0) << " " << sliced_body[2][0][0].M(0) << std::endl;
                break;
            }
            case ADC:{
                ADC_Readout(SeqList[i].t, SeqList[i].line_index, SeqList[i].sample_index, SeqList[i].Gx, SeqList[i].Gy);
                // std::cout << " " <<  body[1][0][0].M(0) << " " << body[1][0][0].M(1) << " " << body[1][0][0].M(2) << std::endl;
                
                break;
            }
        }
    }

}