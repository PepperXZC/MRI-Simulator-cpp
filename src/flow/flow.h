#ifndef FLOW_H
#define FLOW_H
#include "../SeqLoader/SeqLoader.h"
#include "../encoding/encoding.h"
#include "../info/info.h"
#include "omp.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

// Flow 要完成的内容: 将bssfp, molli的yaml序列结合到一起，组成一个新的序列
// 将血池流动关于时间的变化考虑在内，如果流动，就将前两个结点切开，中间插入表示流动的结点
// 做一个 data list 将 8 个矩阵放在一个 vector 里保存 打印到txt文件
// 每张图片对应一个文件夹+2个txt
using std::ofstream;
using std::string;
using std::vector;
using YAML::Node;

typedef struct {
  operation op_node;
  int index;
} generate_node;

class flow_sequence : public molli_sequence {
private:
  // pool_info pool_args;
  double uniform_time;
  double rest_time;
  double each_time_flow;
  bool ro_condition;

public:
  vector<operation> flow_molli_seq;
  flow_sequence(const molli_sequence &msq, pool_info &pi);
  ~flow_sequence();
  void flow_seq_modify();
  void print_flowseq();
  void flow(pool &pl);
};

// flow consideration is combined with molli
class flow_experiment : public Simulator {
private:
  // vector<generate_node> generate_list;
  vector<operation> generate_list;
  int new_proton_index = 1;
  bool reading;
  bool slice_selected = 0;

public:
  flow_experiment(const vector<operation> &seq, pool &pl, const Simulator &sm);
  ~flow_experiment();
  void load_flow_sequence(const vector<operation> &seq, pool &pl);
  void new_proton_generate(pool &pl, int N);
  void save_mat();
  vector<Mat> real_data_list, img_data_list;
};

#endif