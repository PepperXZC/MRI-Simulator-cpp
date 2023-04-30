#include "SeqLoader.h"

map<string, operation_type> operation_type_map = {
    {"NONE", NONE},         {"PULSE", PULSE},     {"ADC", ADC},
    {"ENCODING", ENCODING}, {"READOUT", READOUT}, {"FLOW", FLOW}};

operation::operation() {
  this->FA = 0;
  this->Gx = 0;
  this->Gy = 0;
  this->thickness = 0;
  this->t = 0;
  this->type = NONE;
  this->readout_seq = 0;
}

operation::operation(const operation &op) {
  this->FA = op.FA;
  this->Gx = op.Gx;
  this->Gy = op.Gy;
  this->thickness = op.thickness;
  this->t = op.t;
  this->type = op.type;
  this->line_index = op.line_index;
  this->sample_index = op.sample_index;
  this->readout_seq = op.readout_seq;
}

SeqLoader::SeqLoader(std::string path) {
  this->sequence_path = path;
  this->seq_load();
}

SeqLoader::~SeqLoader() {}

void SeqLoader::seq_load() {
  this->seq_node = YAML::LoadFile(sequence_path);
  for (auto temp : this->seq_node) {
    operation_type type = operation_type_map[temp["type"].as<string>()];
    operation now_info;
    now_info.type = type;
    now_info.t = temp["t"].as<double>();
    now_info.readout_seq = 1;
    switch (type) {
    case NONE:
      break;
    case PULSE: {
      now_info.FA = temp["FA"].as<double>();
      now_info.thickness = temp["slice_thickness"].as<double>();
      break;
    }
    case ENCODING: {
      now_info.Gy = temp["Gy"].as<double>();
      now_info.Gx = temp["Gx"].as<double>();
      break;
    }
    case ADC: {
      now_info.Gx = temp["Gx"].as<double>();
      now_info.line_index = temp["line_index"].as<int>();
      now_info.sample_index = temp["sample_index"].as<int>();
      break;
    }
    }
    this->SeqList.push_back(now_info);
  }
}

template <typename T> vector<T> &operator+(vector<T> &v1, vector<T> &v2) {
  v1.insert(v1.end(), v2.begin(), v2.end());
  return v1;
}

map<string, operation_type> molli_type_map = {
    {"NONE", NONE}, {"PULSE", PULSE}, {"READOUT", READOUT}};

molli_sequence::molli_sequence(const SeqLoader &ro, const string mpath)
    : SeqLoader(ro) {
  molli_path = mpath;
  molli_read();
  // flow = f;
}

molli_sequence::~molli_sequence() {}

void molli_sequence::molli_read() {
  this->molli_node = YAML::LoadFile(molli_path);
  for (auto moll : this->molli_node) {
    operation_type molli_type = molli_type_map[moll["type"].as<string>()];
    operation now_info;
    now_info.type = molli_type;
    if (molli_type != READOUT)
      now_info.t = moll["t"].as<double>();
    // if (now_info.t != 0 && flow == 1) check_flow_time(now_info.t);
    switch (molli_type) {
    case NONE:
      break;
    case PULSE: {
      now_info.FA = moll["FA"].as<double>();
      now_info.thickness = moll["slice_thickness"].as<double>();
      break;
    }
    case READOUT: {
      // this->molli_list = this->molli_list + SeqList;
      now_info.readout_seq = 1;
      break;
    }
    }
    this->molli_list.push_back(now_info);
    if (now_info.readout_seq == 1) {
      this->molli_list = this->molli_list + SeqList;
      operation end_readout;
      end_readout.type = READOUT;
      end_readout.readout_seq = 0;
      this->molli_list.push_back(end_readout);
    }
  }
}