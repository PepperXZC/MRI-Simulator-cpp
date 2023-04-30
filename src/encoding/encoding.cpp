#include "encoding.h"

void Simulator::mat_initialize(int n_pe, int n_read) {
  real_data.resize(n_pe, n_read);
  img_data.resize(n_pe, n_read);
  test_Mat.resize(n_pe, n_read);
}

Simulator::Simulator(pool &sample) {
  // this->pool_ptr = &sample;
  mat_initialize(sample.pool_args.N_pe, sample.pool_args.N_read);
  slice_lower = 0;
  slice_upper = sample.z_length - 1;
}
Simulator::~Simulator() {}

void Simulator::slice_select(pool &sample, double z0, double thickness) {
  int int_z0 = z0 / sample.pool_args.delta;
  int int_tikn = thickness / sample.pool_args.delta;
  slice_lower = int_z0 - int_tikn / 2;
  slice_upper = int_z0 + int_tikn / 2;
}

void Simulator::RF_pulse(pool &sample, double fa) {
  // default x_rot
  Matrix3d fa_mat = Rx(fa * M_PI / 180);
#pragma omp parallel for
  for (int i = 0; i < sample.x_length; i++)
    // #pragma omp parallel for
    for (int j = 0; j < sample.y_length; j++)
      // #pragma omp parallel for
      for (int k = 0; k < sample.z_length; k++) {
        for (int u = 0; u < sample.body[i][j][k].protons.size(); u++) {
          proton_info new_proton(sample.body[i][j][k].flow,
                                 fa_mat * sample.body[i][j][k].protons[u].M);
          sample.body[i][j][k].protons[u] = new_proton;
        }
        sample.body[i][j][k].proton_sum();
      }
  // std::cout << "hi" << std::endl;
}

void Simulator::None_operation(pool &sample, double t) {
  encoding(sample, t, 0, 0, 0);
}

void Simulator::encoding(pool &sample, double t, double Gx = 0, double Gy = 0,
                         double Gz = 0) {
#pragma omp parallel for
  for (int i = 0; i < sample.x_length; i++) {
    for (int j = 0; j < sample.y_length; j++) {
      for (int k = slice_lower; k <= slice_upper; k++) {
        Vector3d pos = sample.body[i][j][k].position;
        double gradient = g * (Gx * pos(0) + Gy * pos(1) + Gz * pos(2));

        Arg_Mats args = freeprecess(t, sample.body[i][j][k].T1,
                                    sample.body[i][j][k].T2, gradient);
        // #pragma omp parallel for
        // if (i == 5 && k == slice_lower) {
        //   std::cout << "position: " << pos(0) << " " << pos(1) << " " <<
        //   pos(2)
        //             << " Gx: " << Gx << " Gy: " << Gy
        //             << " Gradient: " << gradient << std::endl;
        //   std::cout << "proton: " << sample.body[i][j][k].protons[0].M(0) <<
        //   " "
        //             << sample.body[i][j][k].protons[0].M(1) << std::endl;
        //   std::cout << "T: " << sample.body[i][j][k].T1 << " "
        //             << sample.body[i][j][k].T2 << std::endl;
        //   std::cout << "proton A+B: "
        //             << (args.A * sample.body[i][j][k].protons[0].M +
        //             args.B)(0)
        //             << " "
        //             << (args.A * sample.body[i][j][k].protons[0].M +
        //             args.B)(1)
        //             << " "
        //             << (args.A * sample.body[i][j][k].protons[0].M +
        //             args.B)(2)
        //             << std::endl;
        //   std::cout << "voxel: " << sample.body[i][j][k].M(2) << std::endl;
        // }
        // if (i == 32 && k == slice_lower) {
        //   std::cout << "position: " << pos(0) << " " << pos(1) << " " <<
        //   pos(2)
        //             << " Gx: " << Gx << " Gy: " << Gy
        //             << " Gradient: " << gradient << std::endl;
        //   std::cout << "proton: " << sample.body[i][j][k].protons[0].M(0) <<
        //   " "
        //             << sample.body[i][j][k].protons[0].M(1) << std::endl;
        //   std::cout << "T: " << sample.body[i][j][k].T1 << " "
        //             << sample.body[i][j][k].T2 << std::endl;
        //   std::cout << "proton A+B: "
        //             << (args.A * sample.body[i][j][k].protons[0].M +
        //             args.B)(0)
        //             << " "
        //             << (args.A * sample.body[i][j][k].protons[0].M +
        //             args.B)(1)
        //             << " "
        //             << (args.A * sample.body[i][j][k].protons[0].M +
        //             args.B)(2)
        //             << std::endl;
        //   std::cout << "voxel: " << sample.body[i][j][k].M(2) << std::endl;
        // }
        for (int u = 0; u < sample.body[i][j][k].protons.size(); u++) {
          // auto M = args.A * sample.body[i][j][k].protons[u].M + args.B;
          // std::cout << sample.body[i][j][k].protons[u].M << std::endl;
          proton_info new_proton(sample.body[i][j][k].flow,
                                 args.A * sample.body[i][j][k].protons[u].M +
                                     args.B);
          sample.body[i][j][k].protons[u] = new_proton;
        }
        sample.body[i][j][k].proton_sum();
        // if (i == 5 && k == slice_lower) {
        //   std::cout << "voxel A+B: " << sample.body[i][j][k].M(1) <<
        //   std::endl;
        // }
        // if (i == 30) {
        //   std::cout << sample.body[i][0][0].M(0) << " " <<
        //   sample.body[i][0][0].M(1)
        //             << std::endl;
        //   std::cout << sample.body[i][1][0].M(0) << " " <<
        //   sample.body[i][1][0].M(1)
        //             << std::endl;
        //   std::cout << sample.body[i][2][0].M(0) << " " <<
        //   sample.body[i][2][0].M(1)
        //             << std::endl;
        // }
      }
      // std::cout << sample.body[i][j][slice_lower].M(0) << " "
      //           << sample.body[i][j][slice_lower].M(1) << std::endl;
    }
  }
  // std::cout << body[0][0][32].M(1) << " " << body[0][1][32].M(1) << " " <<
  // body[0][2][32].M(1) << std::endl; std::cout << body[1][0][32].M(1) << " "
  // << body[1][1][32].M(1) << " " << body[1][2][32].M(1) <<  std::endl;
  // std::cout << body[2][0][32].M(1) << " " << body[2][1][32].M(1) << " " <<
  // body[2][2][32].M(1) <<  std::endl;
  // std::cout << "hi" << std::endl;
}

void Simulator::ADC_Readout(pool &sample, double t, int line_index,
                            int samsamplee_index, double Gx = 0, double Gy = 0,
                            double Gz = 0) {
  // usually readout after delta_t

  // default: center - started
  double sum_test_data = 0;
  double real_sum = 0, img_sum = 0;
#pragma omp parallel for reduction(+ : real_sum, img_sum, sum_test_data)
  // #pragma omp parallel for reduction(+:sum_test_data)
  for (int i = 0; i < sample.x_length; i++) {
    // #pragma omp parallel for
    for (int j = 0; j < sample.y_length; j++) {
      for (int k = slice_lower; k <= slice_upper; k++) {
        Vector3d pos = sample.body[i][j][k].position;
        double gradient = g * (Gx * pos(0) + Gy * pos(1) + Gz * pos(2));
        Arg_Mats args = freeprecess(t, sample.body[i][j][k].T1,
                                    sample.body[i][j][k].T2, gradient);
        // #pragma omp parallel for
        for (int u = 0; u < sample.body[i][j][k].protons.size(); u++) {
          // std::cout << body[i][j][k].protons[u].M << std::endl;
          proton_info new_proton(sample.body[i][j][k].flow,
                                 args.A * sample.body[i][j][k].protons[u].M +
                                     args.B);

          sample.body[i][j][k].protons[u] = new_proton;
          // std::cout << args.A * body[i][j][k].protons[u].M + args.B <<
          // std::endl;
        }

        sample.body[i][j][k].proton_sum();
        // std::cout << " " <<  body[i][j][k].M(0) << " " << body[i][j][k].M(1)
        // << " " << body[i][j][k].M(2) << std::endl;
        proton_info sum(sample.body[i][j][k].flow, sample.body[i][j][k].M);
        complex<double> Mxy = sum.Mxy;
        // real_sum += Mxy.real() ; img_sum += Mxy.imag();
        real_sum += sample.body[i][j][k].M(0);
        img_sum += sample.body[i][j][k].M(1);
        // if (samsamplee_index % 2 == 0){
        //     real_sum += ( - Mxy.real() ) ; img_sum += ( - Mxy.imag() );
        // }
        // else{
        //     real_sum  += Mxy.real() ; img_sum += Mxy.imag();
        // }
        // only sum the mean of all protons in each voxel
        // sum_test_data = sum_test_data + sum.amsampleitude;
      }
    }
  }

  real_data(line_index, samsamplee_index) = real_sum;
  img_data(line_index, samsamplee_index) = img_sum;
  //   std::cout << "hi" << std::endl;
}

void Simulator::load_seqence(const SeqLoader &Seq, pool &sample) {
  // for (int i = 0; i < tqdm::range(SeqList.size()); i++)
  double time = 0;
  for (int i = 0; i < Seq.SeqList.size(); i++) {
    // bar.progress(i, Seq.SeqList.size());
    switch (Seq.SeqList[i].type) {

    case NONE: {
      None_operation(sample, Seq.SeqList[i].t);
      break;
    }
    case PULSE: {
      RF_pulse(sample, Seq.SeqList[i].FA);
      // std::cout << "step:" << i + 1 << " " << SeqList[i].type << std::endl;
      std::cout << "PULSE!" << std::endl;
      if (slice_selected == 0) {
        slice_select(sample, sample.pool_args.z0, Seq.SeqList[i].thickness);
        slice_selected = 1;
      }
      // std::cout << slice_lower << " " << slice_upper << std::endl;
      // std::cout << " " <<  body[0][0][0].M(0) << " " << body[0][0][0].M(0) <<
      // " " << body[0][0][0].M(0) << std::endl; std::cout << " " <<
      // body[0][1][0].M(0) << " " << body[0][1][0].M(1) << " " <<
      // body[0][1][0].M(2) << std::endl; std::cout << " " <<  body[1][0][0].M
      // << std::endl; std::cout << " " <<  body[1][1][0].M(0) << " " <<
      // body[1][1][0].M(1) << " " << body[1][1][0].M(2) << std::endl;
      break;
    }
    case ENCODING: {
      std::cout << time << std::endl;
      encoding(sample, Seq.SeqList[i].t, Seq.SeqList[i].Gx, Seq.SeqList[i].Gy,
               0);
      break;
    }
    case ADC: {
      time += Seq.SeqList[i].t;
      ADC_Readout(sample, Seq.SeqList[i].t, Seq.SeqList[i].line_index,
                  Seq.SeqList[i].sample_index, Seq.SeqList[i].Gx,
                  Seq.SeqList[i].Gy);
      break;
    }
    }
  }
  // bar.finish();
}