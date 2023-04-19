#include "SeqLoader.h"

map<string, operation_type> operation_type_map = 
{
    {"NONE", NONE},
    {"PULSE", PULSE},
    {"ADC", ADC},
    {"ENCODING", ENCODING}
};

operation::operation()
{
    this->FA = 0;
    this->Gx = 0;
    this->Gy = 0;
    this->thickness = 0;
    this->t = 0;
    this->type = NONE;
}

SeqLoader::SeqLoader(std::string path)
{
    this->path = path;
    this->seq_load();
}

SeqLoader::~SeqLoader()
{
}


void SeqLoader::seq_load()
{
    this->seq_node = YAML::LoadFile(path);
    for (auto temp : this->seq_node)
    {
        operation_type type = operation_type_map[temp["type"].as<string>()];
        operation now_info;
        now_info.type = type;
        now_info.t = temp["t"].as<double>();
        switch (type)
        {
            case NONE:
                break;
            case PULSE:{
                now_info.FA = temp["FA"].as<double>();
                now_info.thickness = temp["slice_thickness"].as<double>();
                break;
            }
            case ENCODING:{
                now_info.Gy = temp["Gy"].as<double>();
                now_info.Gx = temp["Gx"].as<double>();
                break;
            }
            case ADC:{
                now_info.Gx = temp["Gx"].as<double>();
                now_info.line_index = temp["line_index"].as<int>();
                now_info.sample_index = temp["sample_index"].as<int>();
                break;
            }
        }
        this->SeqList.push_back(now_info);
    }
}