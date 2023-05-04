import main_info
import yaml
import numpy as np
import copy
test_info = main_info.info(thickness=0.5)


def sign(num):
    if num > 0:
        return 1
    else:
        return -1


def get_sequence_info(info, prep_num):
    fa_sequence, TR_sequence = np.zeros(
        info.N_pe + prep_num), (np.ones(info.N_pe + prep_num) * info.TR)
    fa_sequence[0], TR_sequence[0] = info.fa / 2, info.TR / 2
    for i in range(1, info.N_pe + + prep_num):
        fa_sequence[i] = sign((i+1) % 2) * info.fa
    return fa_sequence, TR_sequence


def TR_node(fa_sequence, TR_sequence, i, info):
    fa = fa_sequence.tolist()[i]
    t = TR_sequence.tolist()[i]
    init_pulse = {
        "type": "PULSE",
        "FA": fa,
        "t": 0,
        "slice_thickness": info.thickness  # 选层
    }
    res = [init_pulse]
    if i == 0:  # 只给了一个 prep
        free = {
            "type": "NONE",
            "t": t,
        }
        res.append(free)
        return res
    else:
        center = info.N_pe // 2
        num = center + sign((i) % 2) * ((i) // 2)
        G_diff = info.Gyp - num * info.Gyi
        Gx = info.Gx
        phase_encode = {
            "type": "ENCODING",
            "Gy": G_diff,
            "Gx": - Gx,
            "t": info.tau_y
        }

        rewind = {
            "type": "ENCODING",
            "Gy": - G_diff,
            "Gx": - Gx,
            "t": info.tau_y
        }
        res.append(copy.deepcopy(phase_encode))
        for i in range(info.N_read):
            readout = {
                "type": "ADC",
                "Gx": Gx,
                "t": info.delta_t,
                "line_index": num,
                "sample_index": i
            }
            res.append(copy.deepcopy(readout))
        res.append(copy.deepcopy(rewind))
        return res


def molli_generate(info: main_info.info):
    ms_per_beat = 60 * 1e3 / info.HR
    TR_time = info.TR * (info.N_pe + 0.5 + info.prep_num - 1)  # 读取 整整一张图 的时间
    before_lines = info.prep_num - 1  # 第一根线就是中心线
    center_line_time = info.TR * (before_lines + 0.5) + info.TE
    TI_5_before,  TI_3_before = info.TI_5 - \
        center_line_time, info.TI_3 - center_line_time
    interval = ms_per_beat - TR_time
    before_time = [TI_5_before, TI_3_before]
    inversion_interval = 3 * ms_per_beat

    res = []
    readout_time = []
    inversion_pulse = {
        "type": "PULSE",
        "FA": 180,
        "t": 0,
        "slice_thickness": info.thickness  # 选层
    }
    readout_node = {
        "type": "READOUT"
    }
    interval_node = {
        "type": "NONE",
        "t": interval,
    }
    inversion_interval_node = {
        "type": "NONE",
        "t": inversion_interval,
    }
    for i in range(len(before_time)):
        res.append(copy.deepcopy(inversion_pulse))
        time = 0
        node = {
            "type": "NONE",
            "t": before_time[i],
        }
        time += before_time[i]
        res.append(copy.deepcopy(node))
        if i == 0:  # 表示 TI_5
            for j in range(5):
                time += center_line_time
                readout_time.append(time)
                res.append(copy.deepcopy(readout_node))
                time += (TR_time - center_line_time)
                res.append(copy.deepcopy(interval_node))
                time += interval
            time += inversion_interval
            res.append(copy.deepcopy(inversion_interval_node))
        elif i == 1:
            for j in range(3):
                time += center_line_time
                readout_time.append(time)
                res.append(copy.deepcopy(readout_node))
                time += (TR_time - center_line_time)
                res.append(copy.deepcopy(interval_node))
                time += interval
    return res, readout_time


def bssfp_generate(info):
    prep_num = 1
    fa_sequence, TR_sequence = get_sequence_info(info, prep_num)

    seq_list = []
    for i in range(info.N_pe + 1):
        seq_list += TR_node(fa_sequence, TR_sequence, i, info)
    return seq_list


def yaml_generate(info: main_info.info):
    # bssfp
    # seq_list = bssfp_generate(info)
    seq_list, readout_time = molli_generate(info)
    print(readout_time)
    # with open("sequence/molliseq.yaml", "w") as f:
    # with open("sequence/bSSFP_TR2o8.yaml", "w") as f:
    #     yaml.safe_dump(seq_list, f)


yaml_generate(test_info)
# print(sign((0 - 1) % 2) * ((0 + 1) // 2))
