#include "flow.h"
#include <new>
#include <sstream>
#include <vector>

flow_sequence::flow_sequence(const molli_sequence &msq, pool_info &pi)
    : molli_sequence(msq) {
  uniform_time = 0;
  each_time_flow = pi.delta / (pi.flow_speed * 1e-3);
  ro_condition = 0;
  flow_seq_modify();
}

flow_sequence::~flow_sequence() {}

void flow_sequence::flow_seq_modify() {
  for (operation op_node : molli_list) {
    if (op_node.type == READOUT)
      ro_condition = op_node.readout_seq;
    if (op_node.t != 0) {
      if (ro_condition == 0) {
        operation copy_node(op_node);
        flow_molli_seq.push_back(copy_node);
        uniform_time = fmod((op_node.t + uniform_time), each_time_flow);
        continue;
      }

      double new_ut = 0;
      double rest_time = 0;
      int num_flow = (op_node.t + uniform_time) / each_time_flow;
      bool flag = 0;
      double rt = fmod((op_node.t + uniform_time), each_time_flow);
      if (rt == (op_node.t + uniform_time)) {
        rest_time = op_node.t;
        uniform_time += op_node.t;
      } else {
        rest_time = rt;
        new_ut = rest_time;
        flag = 1;
      }

      for (int n = 0; n < num_flow; n++) {
        operation flow_node;
        flow_node.type = FLOW;
        double t;
        if (n == 0)
          t = each_time_flow - uniform_time;
        else
          t = each_time_flow;
        operation forward_half_node(op_node);
        forward_half_node.t = t;
        flow_molli_seq.push_back(forward_half_node);
        flow_molli_seq.push_back(flow_node);
      }
      operation end_half_node(op_node);
      end_half_node.t = rest_time;
      flow_molli_seq.push_back(end_half_node);

      if (flag == 1)
        uniform_time = new_ut;
    } else {
      operation copy_node(op_node);
      flow_molli_seq.push_back(copy_node);
    }
  }
}

void flow_sequence::print_flowseq() {
  // for (operation flow_node : flow_molli_seq) {
  double time = 0;
  bool reading = 0;
  for (int i = 0; i < flow_molli_seq.size(); i++) {
    if (flow_molli_seq[i].type == READOUT && flow_molli_seq[i].readout_seq == 1)
      reading = 1;
    else if (flow_molli_seq[i].type == READOUT &&
             flow_molli_seq[i].readout_seq == 0)
      reading = 0;
    if (reading == 1)
      time += flow_molli_seq[i].t;
    // std::cout << flow_molli_seq[i].type << " " << flow_molli_seq[i].t << "
    // "
    //           << std::endl;
    if (flow_molli_seq[i].type == FLOW) {
      std::cout << time << " " << (time == 5) << std::endl;
      time = 0;
    }
  }
}

flow_experiment::flow_experiment(const vector<operation> &seq, pool &pl,
                                 const Simulator &sm)
    : Simulator(sm) {
  // print_flowseq();
  load_flow_sequence(seq, pl);
}

flow_experiment::~flow_experiment() {}

void flow_experiment::load_flow_sequence(const vector<operation> &seq,
                                         pool &pl) {
  reading = 0;
  tqdm bar;
  int num_Mz = 0;
  for (int i = 0; i < seq.size(); i++) {
    bar.progress(i, seq.size());
    operation_type type = seq[i].type;
    switch (type) {
    case NONE:
      None_operation(pl, seq[i].t);
      break;
    case PULSE: {
      RF_pulse(pl, seq[i].FA);
      if (slice_selected == 0) {
        slice_select(pl, pl.pool_args.z0, seq[i].thickness);
        slice_selected = 1;
      }
      break;
    }
    case READOUT: {
      if (seq[i].readout_seq == 1)
        reading = 1;
      else {
        reading = 0;
        // real_data_list.push_back(real_data);
        // img_data_list.push_back(img_data);
        real_data_list.push_back(real_data);
        img_data_list.push_back(img_data);
        ofstream fout;
        fout.open("real_mat.txt", std::ios::trunc);
        fout << real_data << std::endl;
        fout.close();
        fout.open("img_mat.txt", std::ios::trunc);
        fout << img_data << std::endl;
        fout.close();
        std::cout << "Reaout!" << std::endl;
        Mat temp_Mz;
        temp_Mz.resize(pl.pool_args.N_pe, pl.pool_args.N_read);
        for (int i = 0; i < Mz_list.size(); i++) {
          temp_Mz += Mz_list[i];
        }
        temp_Mz /= num_Mz;
        fout.open("mat_data.txt", std::ios::trunc);
        fout << temp_Mz << std::endl;
        fout.close();
        // Mz_data_list.push_back(Mz_list);
        Mz_data_list.push_back(temp_Mz);
        Mz_list.clear();
      }
      break;
    }

    case ENCODING: {
      encoding(pl, seq[i].t, seq[i].Gx, seq[i].Gy, 0);
      break;
    }

    case ADC: {
      ADC_Readout(pl, seq[i].t, seq[i].line_index, seq[i].sample_index,
                  seq[i].Gx, seq[i].Gy, 0);
      num_Mz += 1;
      break;
    }
    case FLOW: {
      for (int i = 0; i < pl.pool_args.num_vassels; i++) {
        pl.pool_roll(i);
        new_proton_generate(i, pl, new_proton_index);
      }
      new_proton_index += 1;
      break;
    }
    }
    // if (flow_molli_seq[i].t != 0)
    generate_list.push_back(seq[i]);
  }
  bar.finish();
}

void flow_experiment::new_proton_generate(int vassel_index,
                                          const pool &pl_example, int N = 0) {
  double y_max = pl_example.pool_args.fov / 2;
  double y0 = 0;
  if (N != 0)
    y0 = y_max + N * pl_example.pool_args.delta;
  int pool_width =
      pl_example.upper_list[vassel_index] - pl_example.lower_list[vassel_index];
  pool temp_pool(pl_example.pool_args, pool_width * pl_example.pool_args.delta,
                 pl_example.pool_args.delta, pl_example.pool_args.fov);

  // Voxel **new_voxels = new Voxel *[pool_width];
  // 这里是较为标准化的 pool init 操作
  // #pragma omp parallel for
  temp_pool.vassel_init(pl_example.pool_args.T_vassel[vassel_index][0],
                        pl_example.pool_args.T_vassel[vassel_index][1],
                        pl_example.lower_list[vassel_index], y0);
  for (int i = 0; i < generate_list.size(); i++) {
    operation_type type = generate_list[i].type;
    switch (type) {
    case NONE:
      None_operation(temp_pool, generate_list[i].t);
      break;
    case PULSE:
      RF_pulse(temp_pool, generate_list[i].FA);
      // std::cout << "                " <<
      // temp_pool.body[5][0][slice_lower].M(1)
      //           << " " << std::endl;

      break;
    case READOUT:
      break;
    case ENCODING: {
      encoding(temp_pool, generate_list[i].t, generate_list[i].Gx,
               generate_list[i].Gy, 0);
      break;
    }
    case ADC:
      encoding(temp_pool, generate_list[i].t, generate_list[i].Gx,
               generate_list[i].Gy, 0);
      break;
    case FLOW: {
#pragma omp parallel for
      for (int p = 0; p < temp_pool.x_length; p++) {
        for (int j = 0; j < temp_pool.y_length; j++) {
          for (int k = 0; k < temp_pool.z_length; k++) {
            // temp_pool.y_length == 1
            // Vector3d pos = temp_pool.body[p][j][k].position;
            Vector3d pos = {temp_pool.body[p][j][k].position(0),
                            temp_pool.body[p][j][k].position(1) -
                                pl_example.pool_args.delta,
                            temp_pool.body[p][j][k].position(2)};
            temp_pool.body[p][j][k].initialize(
                pl_example.pool_args.T_vassel[vassel_index][0],
                pl_example.pool_args.T_vassel[vassel_index][1], 1, pos,
                temp_pool.body[p][j][k].M, 1);
          }
        }
      }
      break;
    }
    }
  }

#pragma omp parallel for
  for (int i = pl_example.lower_list[vassel_index];
       i < pl_example.upper_list[vassel_index]; i++)
    for (int k = 0; k < pl_example.z_length; k++)
      // pl.body[i][0][k].initialize(temp_pool.body[i - lower][0][k]);
      pl_example.body[i][0][k].initialize(
          pl_example.pool_args.T_vassel[vassel_index][0],
          pl_example.pool_args.T_vassel[vassel_index][1], 1,
          pl_example.body[i][0][k].position,
          temp_pool.body[i - pl_example.lower_list[vassel_index]][0][k].M, 1);
}

void flow_experiment::save_mat(const string &path) {
  for (int i = 0; i < real_data_list.size(); i++) {
    // std::cout << test_simulator.test_Mat << std::endl;
    std::ostringstream oss;
    oss << i;
    string num_str = oss.str();
    string r_txt = "real_mat_";
    r_txt = path + r_txt + num_str + ".txt";
    string i_txt = "img_mat_";
    i_txt = path + i_txt + num_str + ".txt";
    string z_txt = "Mz_mat_";
    z_txt = path + z_txt + num_str + ".txt";
    ofstream fout;
    std::cout << "path: " << r_txt << " " << i_txt << std::endl;
    fout.open(r_txt, std::ios::trunc);
    fout << real_data_list[i] << std::endl;
    fout.close();
    fout.open(i_txt, std::ios::trunc);
    fout << img_data_list[i] << std::endl;
    fout.close();
    fout.open(z_txt, std::ios::trunc);
    fout << Mz_data_list[i] << std::endl;
    fout.close();
  }
}
