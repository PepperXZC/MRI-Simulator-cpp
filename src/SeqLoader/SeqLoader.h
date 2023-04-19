#ifndef SEQLOADER_H
#define SEQLOADER_H
#include <string>
#include <map>
#include <iostream>
#include <yaml-cpp/yaml.h>

using std::map;
using std::string;
using std::vector;
using YAML::Node;

enum operation_type
{
    NONE,
    PULSE,
    ADC,
    ENCODING
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
    operation();
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
    std::string path;
};
#endif