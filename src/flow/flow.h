#ifndef FLOW_H
#define FLOW_H
#include <iostream>
#include "../info/info.h"
#include "../SeqLoader/SeqLoader.h"
#include "../info/info.h"
#include "../encoding/encoding.h"

// Flow 要完成的内容: 将bssfp, molli的yaml序列结合到一起，组成一个新的序列
// 将血池流动关于时间的变化考虑在内，如果流动，就将前两个结点切开，中间插入表示流动的结点
// 做一个 data list 将 8 个矩阵放在一个 vector 里保存 打印到txt文件 每张图片对应一个文件夹+2个txt
using std::string;
using std::vector;
using YAML::Node;

class flow_sequence: public molli_sequence
{
private:
    pool_info pool_args;
    double uniform_time;
    double rest_time;
    double each_time_flow;
    bool ro_condition;
public:
    vector<operation> flow_molli_seq;
    flow_sequence(const molli_sequence& msq, pool_info& pi);
    ~flow_sequence();
    void flow_seq_modify();
    void print_flowseq();
    void flow(pool& pl);
};
// flow consideration is combined with molli

class flow_experiment : public flow_sequence, pool, Simulator
{
private:
public:
    flow_experiment(const flow_sequence& fs, pool& pl, const Simulator& sm);
    ~flow_experiment();
    void load_flow_sequence();
};

#endif