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
      break;
    }
    case FLOW: {
      pl.pool_roll();
      new_proton_generate(pl, new_proton_index);
      new_proton_index += 1;
      break;
    }
    }
    // if (flow_molli_seq[i].t != 0)
    generate_list.push_back(seq[i]);
  }
  bar.finish();
}

void flow_experiment::new_proton_generate(pool &pl, int N = 0) {
  double y_max = pl.pool_args.fov / 2;
  double y0 = 0;
  if (N != 0)
    y0 = y_max + N * pl.pool_args.delta;

  int pool_width = pl.upper - pl.lower;
  pool temp_pool(pl.pool_args, pool_width * pl.pool_args.delta,
                 pl.pool_args.delta, pl.pool_args.fov);

  // Voxel **new_voxels = new Voxel *[pool_width];
  // 这里是较为标准化的 pool init 操作
  temp_pool.body = new Voxel **[temp_pool.x_length];
#pragma omp parallel for
  for (int i = 0; i < temp_pool.x_length; i++) {
    temp_pool.body[i] = new Voxel *[temp_pool.y_length];
    for (int j = 0; j < temp_pool.y_length; j++) {
      temp_pool.body[i][j] = new Voxel[temp_pool.z_length];
      for (int k = 0; k < temp_pool.z_length; k++) {
        Vector3d M_init;
        M_init << 0, 0, 1;
        Vector3d pos;
        pos << -pl.pool_args.fov / 2 + (i + pl.lower) * pl.pool_args.delta, y0,
            k * pl.pool_args.delta;
        temp_pool.body[i][j][k].initialize(pl.pool_args.T_vassel, 1, pos,
                                           M_init, 1, 0);
      }
    }
  }
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
            Vector3d pos;
            pos << temp_pool.body[p][j][k].position(0),
                temp_pool.body[p][j][k].position(1) - pl.pool_args.delta,
                temp_pool.body[p][j][k].position(2);
            temp_pool.body[p][j][k].initialize(pl.pool_args.T_vassel, 1, pos,
                                               temp_pool.body[p][j][k].M, 1, 0);
            // std::cout << temp_pool.body[p][j][k].M(0) << " "
            //           << temp_pool.body[p][j][k].M(1) << " "
            //           << temp_pool.body[p][j][k].M(2) << std::endl;
          }
        }
      }
      // std::cout << "New, after flow: "
      //           << temp_pool.body[5][0][slice_lower].position(0) << " "
      //           << temp_pool.body[5][0][slice_lower].position(1) <<
      //           std::endl;
      break;
    }
      // #pragma omp parallel for
    }
  }
  // std::cout << "generate: "
  //           << temp_pool.body[31 - pl.lower][0][slice_lower].M(1) << " "
  //           << temp_pool.body[32 - pl.lower][0][slice_lower].M(1) << " "
  //           << temp_pool.body[33 - pl.lower][0][slice_lower].M(1) <<
  //           std::endl;
#pragma omp parallel for
  for (int i = pl.lower; i < pl.upper; i++)
    for (int k = 0; k < pl.z_length; k++)
      // pl.body[i][0][k].initialize(temp_pool.body[i - lower][0][k]);
      pl.body[i][0][k].initialize(pl.pool_args.T_vassel, 1,
                                  pl.body[i][0][k].position,
                                  temp_pool.body[i - pl.lower][0][k].M, 1, 0);
  // std::cout << "start of pool: " << pl.body[31][0][0].M(1) << " "
  //           << pl.body[32][0][slice_lower].M(1) << " "
  //           << pl.body[33][0][slice_lower].M(1) << std::endl;
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
    ofstream fout;
    std::cout << "path: " << r_txt << " " << i_txt << std::endl;
    fout.open(r_txt, std::ios::trunc);
    fout << real_data_list[i] << std::endl;
    fout.close();
    fout.open(i_txt, std::ios::trunc);
    fout << img_data_list[i] << std::endl;
    fout.close();
  }
}
