#include "encoding.h"

void Simulator::mat_initialize(int n_pe, int n_read)
{
    real_data.resize(n_pe, n_read);
    img_data.resize(n_pe, n_read);
    test_Mat.resize(n_pe, n_read);
}

Simulator::Simulator(pool &pl)
{
    // this->pool_ptr = &pl;
    mat_initialize(pl.pool_args.N_pe, pl.pool_args.N_read);
}
Simulator::~Simulator()
{
}

void Simulator::slice_select(pool &pl, double z0, double thickness)
{
    int int_z0 = z0 / pl.pool_args.delta; int int_tikn = thickness / pl.pool_args.delta;
    slice_lower = int_z0 - int_tikn / 2; slice_upper = int_z0 + int_tikn / 2;
}

void Simulator::RF_pulse(pool &pl, double fa)
{
    // default x_rot
    Matrix3d fa_mat = Rx(fa * M_PI / 180);
    #pragma omp parallel for
    for (int i = 0; i < pl.pool_length; i++)
        // #pragma omp parallel for
        for (int j = 0; j < pl.pool_length; j++)
            // #pragma omp parallel for
            for (int k = 0; k < pl.pool_length; k++){
                for (int u = 0; u < pl.body[i][j][k].protons.size(); u++){
                    proton_info new_proton(pl.body[i][j][k].flow, fa_mat * pl.body[i][j][k].protons[u].M);
                    pl.body[i][j][k].protons[u] = new_proton;
                } 
                pl.body[i][j][k].proton_sum();
            } 
    std::cout << "hi" << std::endl;         
}

void Simulator::None_operation(pool &pl, double t)
{
    encoding(pl, t, 0, 0, 0);
}

void Simulator::encoding(pool &pl, double t, double Gx = 0, double Gy = 0, double Gz = 0)
{
    #pragma omp parallel for
    for (int i = 0; i < pl.pool_length; i++){
        for (int j = 0; j < pl.pool_length; j++){
            for (int k = slice_lower; k <= slice_upper; k++){
                Vector3d pos = pl.body[i][j][k].position;
                double gradient = g * (Gx * pos(0) + Gy * pos(1) + Gz * pos(2));
                
                Arg_Mats args = freeprecess(t, pl.body[i][j][k].T1, pl.body[i][j][k].T2, gradient);
                // #pragma omp parallel for
                for (int u = 0; u < pl.body[i][j][k].protons.size(); u++){
                    proton_info new_proton(pl.body[i][j][k].flow, args.A * pl.body[i][j][k].protons[u].M + args.B);
                    pl.body[i][j][k].protons[u] = new_proton;
                }
                pl.body[i][j][k].proton_sum();
                // if ( i == 1 && j == 0 && k == 0) {std::cout << Gx  << " " << Gy  << " " << Gz  << " " << pl.body[i][j][k].M(0) << " " << body[i][j][k].M(1) << " " << body[i][j][k].M(2) << std::endl;}
            }
        }
    }
    // std::cout << body[0][0][32].M(1) << " " << body[0][1][32].M(1) << " " << body[0][2][32].M(1) << std::endl;
    // std::cout << body[1][0][32].M(1) << " " << body[1][1][32].M(1) << " " << body[1][2][32].M(1) <<  std::endl;
    // std::cout << body[2][0][32].M(1) << " " << body[2][1][32].M(1) << " " << body[2][2][32].M(1) <<  std::endl;
    std::cout << "hi" << std::endl;
}

void Simulator::ADC_Readout(pool &pl, double t, int line_index, int sample_index, double Gx = 0, double Gy = 0, double Gz = 0)
{
    // usually readout after delta_t

    // default: center - started
    double sum_test_data = 0;
    double real_sum = 0, img_sum = 0;
    #pragma omp parallel for reduction(+:real_sum, img_sum, sum_test_data)
    // #pragma omp parallel for reduction(+:sum_test_data)
    for (int i = 0; i < pl.pool_length; i++){
        // #pragma omp parallel for
        for (int j = 0; j < pl.pool_length; j++){
            for (int k = slice_lower; k <= slice_upper; k++){
                Vector3d pos = pl.body[i][j][k].position;
                double gradient = g * (Gx * pos(0) + Gy * pos(1) + Gz * pos(2));
                Arg_Mats args = freeprecess(t, pl.body[i][j][k].T1, pl.body[i][j][k].T2, gradient);
                // #pragma omp parallel for
                for (int u = 0; u < pl.body[i][j][k].protons.size(); u++){
                    // std::cout << body[i][j][k].protons[u].M << std::endl;
                    proton_info new_proton(pl.body[i][j][k].flow, args.A * pl.body[i][j][k].protons[u].M + args.B);
                    
                    pl.body[i][j][k].protons[u] = new_proton;
                    // std::cout << args.A * body[i][j][k].protons[u].M + args.B << std::endl;
                }
                    
                pl.body[i][j][k].proton_sum();
                // std::cout << " " <<  body[i][j][k].M(0) << " " << body[i][j][k].M(1) << " " << body[i][j][k].M(2) << std::endl;
                proton_info sum(pl.body[i][j][k].flow, pl.body[i][j][k].M);
                complex<double> Mxy = sum.Mxy;
                // real_sum += Mxy.real() ; img_sum += Mxy.imag();
                real_sum += pl.body[i][j][k].M(0); img_sum += pl.body[i][j][k].M(1); 
                // if (sample_index % 2 == 0){
                //     real_sum += ( - Mxy.real() ) ; img_sum += ( - Mxy.imag() );  
                // }
                // else{
                //     real_sum  += Mxy.real() ; img_sum += Mxy.imag();
                // }
                     // only sum the mean of all protons in each voxel
                // sum_test_data = sum_test_data + sum.amplitude;
            }
        }
    }
            
    real_data(line_index, sample_index) = real_sum;
    img_data(line_index, sample_index) = img_sum;
    std::cout << "hi" << std::endl;
}

void Simulator::load_seqence(const SeqLoader &Seq, pool &pl)
{
    // for (int i = 0; i < tqdm::range(SeqList.size()); i++)
    double time = 0;
    for (int i = 0; i < Seq.SeqList.size(); i++)
    {
        switch (Seq.SeqList[i].type)
        {
            
            case NONE: {
                None_operation(pl, Seq.SeqList[i].t);
                break;
            }
            case PULSE:{
                RF_pulse(pl, Seq.SeqList[i].FA);
                // std::cout << "step:" << i + 1 << " " << SeqList[i].type << std::endl;
                std::cout << "PULSE!" << std::endl;
                if (slice_selected == 0){
                    slice_select(pl, pl.pool_args.z0, Seq.SeqList[i].thickness);
                    slice_selected = 1;
                }
                // std::cout << slice_lower << " " << slice_upper << std::endl;
                // std::cout << " " <<  body[0][0][0].M(0) << " " << body[0][0][0].M(0) << " " << body[0][0][0].M(0) << std::endl;
                // std::cout << " " <<  body[0][1][0].M(0) << " " << body[0][1][0].M(1) << " " << body[0][1][0].M(2) << std::endl;
                // std::cout << " " <<  body[1][0][0].M << std::endl;
                // std::cout << " " <<  body[1][1][0].M(0) << " " << body[1][1][0].M(1) << " " << body[1][1][0].M(2) << std::endl;
                break;
            }
            case ENCODING:{
                std::cout << time << std::endl;
                encoding(pl, Seq.SeqList[i].t, Seq.SeqList[i].Gx, Seq.SeqList[i].Gy, 0);
                break;
            }
            case ADC:{
                time += Seq.SeqList[i].t;
                ADC_Readout(pl, Seq.SeqList[i].t, Seq.SeqList[i].line_index, Seq.SeqList[i].sample_index, Seq.SeqList[i].Gx, Seq.SeqList[i].Gy);
                break;
            }
        }
    }

}