#ifndef SEQLOADER_H
#define SEQLOADER_H
#include <string>
#include <map>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include <math.h>

using std::map;
using std::string;
using std::vector;
using YAML::Node;

enum operation_type
{
    NONE,
    PULSE,
    ADC,
    ENCODING,
    READOUT,
    FLOW
};

struct operation
{
    operation_type type;
    double FA;
    double thickness;
    double Gx;          // frequency encoding
    double Gy;          // phase encoding
    // bool readout;       // after this interval t, operate readout process.
    double t;
    int line_index;
    int sample_index;
    bool readout_seq;
    operation();
    operation(const operation &op);
};

class SeqLoader
{
public:
    SeqLoader(std::string path);
    ~SeqLoader();
    vector<operation> SeqList;
private:
    void seq_load();
    Node seq_node;
    std::string sequence_path;
};

class molli_sequence : public SeqLoader
{
private:
    string molli_path;
    Node molli_node;
    // bool flow;
public:
    vector<operation> molli_list;
    void molli_read();
    molli_sequence(const SeqLoader &ro, const string mpath);
    ~molli_sequence();
};
#endif