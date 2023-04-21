#include "encoding.h"

void Simulator::mat_initialize()
{
    real_data.resize(pool_args.N_pe, pool_args.N_read);
    img_data.resize(pool_args.N_pe, pool_args.N_read);
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
    // body = new Voxel **[pool_length];
    int int_z0 = z0 / pool_args.delta; int int_tikn = thickness / pool_args.delta;
    slice_lower = int_z0 - int_tikn / 2; slice_upper = int_z0 + int_tikn / 2;
}

void Simulator::RF_pulse(double fa)
{
    // default x_rot
    Matrix3d fa_mat = Rx(fa * M_PI / 180);
    // #pragma omp parallel
    // #pragma omp for
    // #pragma omp parallel for
    for (int i = 0; i < pool_length; i++)
        // #pragma omp parallel for
        for (int j = 0; j < pool_length; j++)
            // #pragma omp parallel for
            for (int k = 0; k < pool_length; k++){
                for (int u = 0; u < body[i][j][k].protons.size(); u++){
                    proton_info new_proton(body[i][j][k].flow, fa_mat * body[i][j][k].protons[u].M);
                    body[i][j][k].protons[u] = new_proton;
                } 
                body[i][j][k].proton_sum();
            } 
    // std::cout << body[0][0][32].M(1) << " " << body[0][1][32].M(1) << " " << body[0][2][32].M(1) << std::endl;
    // std::cout << body[1][0][32].M(1) << " " << body[1][1][32].M(1) << " " << body[1][2][32].M(1) <<  std::endl;
    // std::cout << body[2][0][32].M(1) << " " << body[2][1][32].M(1) << " " << body[2][2][32].M(1) <<  std::endl;
    std::cout << "hi" << std::endl;         
}

void Simulator::None_operation(double t)
{
    encoding(t, 0, 0, 0);
}

void Simulator::encoding(double t, double Gx = 0, double Gy = 0, double Gz = 0)
{
    // #pragma omp parallel
    // #pragma omp for
    // #pragma omp parallel for
    // std::cout << "hi" << std::endl;
    for (int i = 0; i < pool_length; i++){
        // #pragma omp parallel for
        for (int j = 0; j < pool_length; j++){
            // #pragma omp parallel for
            // std::cout << "hi" << std::endl;
            for (int k = slice_lower; k <= slice_upper; k++){
                Vector3d pos = body[i][j][k].position;
                double gradient = g * (Gx * pos(0) + Gy * pos(1) + Gz * pos(2));
                
                Arg_Mats args = freeprecess(t, body[i][j][k].T1, body[i][j][k].T2, gradient);
                // #pragma omp parallel for
                for (int u = 0; u < body[i][j][k].protons.size(); u++){
                    proton_info new_proton(body[i][j][k].flow, args.A * body[i][j][k].protons[u].M + args.B);
                    // std::cout << body[i][j][k].protons[u].M << std::endl;
                    // std::cout << new_proton.M << std::endl;
                    body[i][j][k].protons[u] = new_proton;
                }
                body[i][j][k].proton_sum();
                if ( i == 1 && j == 0 && k == 0) {std::cout << Gx  << " " << Gy  << " " << Gz  << " " << body[i][j][k].M(0) << " " << body[i][j][k].M(1) << " " << body[i][j][k].M(2) << std::endl;}
            }
        }
    }
    // std::cout << body[0][0][32].M(1) << " " << body[0][1][32].M(1) << " " << body[0][2][32].M(1) << std::endl;
    // std::cout << body[1][0][32].M(1) << " " << body[1][1][32].M(1) << " " << body[1][2][32].M(1) <<  std::endl;
    // std::cout << body[2][0][32].M(1) << " " << body[2][1][32].M(1) << " " << body[2][2][32].M(1) <<  std::endl;
    std::cout << "hi" << std::endl;
}

void Simulator::ADC_Readout(double t, int line_index, int sample_index, double Gx = 0, double Gy = 0, double Gz = 0)
{
    // usually readout after delta_t

    // default: center - started
    double sum_test_data = 0;
    // complex<double> sum_data = {0, 0};
    double real_sum = 0, img_sum = 0;
    // #pragma omp parallel for reduction(+:real_sum, img_sum, sum_test_data)
    // #pragma omp parallel for reduction(+:sum_test_data)
    for (int i = 0; i < pool_length; i++){
        // #pragma omp parallel for
        for (int j = 0; j < pool_length; j++){
            // #pragma omp parallel for
            for (int k = slice_lower; k <= slice_upper; k++){
                Vector3d pos = body[i][j][k].position;
                double gradient = g * (Gx * pos(0) + Gy * pos(1) + Gz * pos(2));
                // if ( i == 1 && j == 0 && k == 0) std::cout << Gx << " " << Gy << " " << gradient << std::endl;
                Arg_Mats args = freeprecess(t, body[i][j][k].T1, body[i][j][k].T2, gradient);
                // if ( i >= 27 && i <= 36)
                //     std::cout << pos(0)  << " " << pos(1)  << " " << pos(2)  << " " << gradient << " " <<  body[i][j][k].T1 << std::endl;
                // std::cout << Gx  << " " << Gy  << " " << Gz  << std::endl;
                // #pragma omp parallel for
                for (int u = 0; u < body[i][j][k].protons.size(); u++){
                    // std::cout << body[i][j][k].protons[u].M << std::endl;
                    proton_info new_proton(body[i][j][k].flow, args.A * body[i][j][k].protons[u].M + args.B);
                    
                    body[i][j][k].protons[u] = new_proton;
                    // std::cout << args.A * body[i][j][k].protons[u].M + args.B << std::endl;
                }
                    
                body[i][j][k].proton_sum();
                // std::cout << " " <<  body[i][j][k].M(0) << " " << body[i][j][k].M(1) << " " << body[i][j][k].M(2) << std::endl;
                // proton_info sum(body[i][j][k].flow, body[i][j][k].M);
                // complex<double> Mxy = sum.Mxy;
                // real_sum += Mxy.real() ; img_sum += Mxy.imag();
                real_sum += body[i][j][k].M(0); img_sum += body[i][j][k].M(1); 
                // if (sample_index % 2 == 0){
                //     real_sum = real_sum + ( - Mxy.real() ) ; img_sum = img_sum + ( - Mxy.imag() );  
                // }
                // else{
                //     real_sum = real_sum + Mxy.real() ; img_sum = img_sum + Mxy.imag();
                // }
                     // only sum the mean of all protons in each voxel
                // sum_test_data = sum_test_data + sum.amplitude;
                // std::cout << Mxy << sum.amplitude << std::endl;
                // if ( i == 1 && j == 0 && k == 0) {std::cout << Gx  << " " << Gy  << " " << Gz  << " " << body[i][j][k].M(0) << " " << body[i][j][k].M(1) << " " << body[i][j][k].M(2) << std::endl;}
            }
        }
    }
            
    real_data(line_index, sample_index) = real_sum;
    img_data(line_index, sample_index) = img_sum;
    // std::cout << body[0][0][32].M(0) << " " << body[63][0][32].M(0) << std::endl;
    // std::cout << body[0][0][32].M(1) << " " << body[63][0][32].M(1) << std::endl;
    // std::cout << body[2][0][32].M(1) << " " << body[2][1][32].M(1) << " " << body[2][2][32].M(1) <<  std::endl;
    std::cout << "hi" << std::endl;
    // std::cout << real_sum << " " << img_sum << std::endl;
    // test_Mat(line_index, sample_index) = sum_test_data;
}

void Simulator::check_readout_record(int readout_index)
{
    // if (readout_index == pool_args.N_read - 1)
    // {
    //     result_list.push_back(data);
    //     mat_initialize();
    // }
}

void Simulator::load_seqence()
{
    // for (int i = 0; i < tqdm::range(SeqList.size()); i++)
    double time = 0;
    for (int i = 0; i < SeqList.size(); i++)
    {
        // std::cout << body[1][0][0].M << std::endl;
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
                // std::cout << slice_lower << " " << slice_upper << std::endl;
                // std::cout << " " <<  body[0][0][0].M(0) << " " << body[0][0][0].M(0) << " " << body[0][0][0].M(0) << std::endl;
                // std::cout << " " <<  body[0][1][0].M(0) << " " << body[0][1][0].M(1) << " " << body[0][1][0].M(2) << std::endl;
                // std::cout << " " <<  body[1][0][0].M << std::endl;
                // std::cout << " " <<  body[1][1][0].M(0) << " " << body[1][1][0].M(1) << " " << body[1][1][0].M(2) << std::endl;
                break;
            }
            case ENCODING:{
                std::cout << time << std::endl;
                encoding(SeqList[i].t, SeqList[i].Gx, SeqList[i].Gy, 0);
                // std::cout << body[0][0][32].M(1) << " " << body[0][1][32].M(1) << std::endl;
                // std::cout << body[1][0][32].M(1) << " " << body[1][1][32].M(1) << std::endl;
                // std::cout << " " <<  body[1][0][32].M(0) << " " << body[1][0][32].M(1) << " " << body[1][0][32].M(2) << std::endl;
                // std::cout << " " <<  body[1][1][32].M(0) << " " << body[1][1][32].M(1) << " " << body[1][1][32].M(2) << std::endl;
                // std::cout << " " <<  body[0][0][0].M(0) << " " << body[1][0][0].M(0) << " " << body[2][0][0].M(0) << std::endl;
                break;
            }
            case ADC:{
                time += SeqList[i].t;
                ADC_Readout(SeqList[i].t, SeqList[i].line_index, SeqList[i].sample_index, SeqList[i].Gx, SeqList[i].Gy);
                // std::cout << " " <<  body[1][0][0].M(0) << " " << body[1][0][0].M(1) << " " << body[1][0][0].M(2) << std::endl;
                if (SeqList[i].sample_index == 31){
                    // std:: cout << "no" << std::endl;
                    // std::cout << " " <<  body[1][0][32].M(0) << " " << body[1][0][32].M(1) << " " << body[1][0][32].M(2) << std::endl;
                    // std::cout << " " <<  body[1][1][32].M(0) << " " << body[1][1][32].M(1) << " " << body[1][1][32].M(2) << std::endl;
                    // std::cout << " " <<  body[0][1][32].M(0) << " " << body[0][1][32].M(1) << " " << body[0][1][32].M(2) << std::endl;
                }
                break;
            }
        }
    }

}