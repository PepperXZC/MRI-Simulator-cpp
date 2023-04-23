#include"flow.h"

flow_sequence::flow_sequence(const molli_sequence& msq, pool_info& pi)
    :molli_sequence(msq)
{
    uniform_time = 0;
    each_time_flow = pi.delta / (pi.flow_speed * 1e-3);
    ro_condition = 0;
    flow_seq_modify();
}

flow_sequence::~flow_sequence(){ }

void flow_sequence::flow_seq_modify()
{
    for (operation op_node : molli_list )
    {
        if (op_node.type == READOUT) ro_condition = op_node.readout_seq;
        if (op_node.t != 0)
        {
            if (ro_condition == 0)
            {
                operation copy_node(op_node);
                flow_molli_seq.push_back(copy_node);
                continue;
            }

            double new_ut = 0; double rest_time = 0; int num_flow = (op_node.t + uniform_time) / each_time_flow;
            bool flag = 0;
            double rt = fmod((op_node.t + uniform_time), each_time_flow);
            if ( rt == (op_node.t + uniform_time)){
                rest_time = op_node.t; 
                uniform_time += op_node.t;
            }
            else{
                rest_time = rt;
                new_ut = rest_time;
                flag = 1;
            }

            for (int n = 0; n < num_flow; n++ )
            {
                operation flow_node; flow_node.type = FLOW; double t;
                if (n == 0) t = each_time_flow - uniform_time; else t = each_time_flow;
                operation forward_half_node(op_node); forward_half_node.t = t;
                flow_molli_seq.push_back(forward_half_node);
                flow_molli_seq.push_back(flow_node);
            }
            operation end_half_node(op_node); end_half_node.t = rest_time;
            flow_molli_seq.push_back(end_half_node);

            if (flag == 1) uniform_time = new_ut;
        }
        else
        {
            operation copy_node(op_node);
            flow_molli_seq.push_back(copy_node);
        }
    }
}

void flow_sequence::print_flowseq()
{
    for (operation flow_node : flow_molli_seq )
    {
        std::cout << flow_node.type << " " << flow_node.t << " ";
        // std::cout << std::endl;
     }
}