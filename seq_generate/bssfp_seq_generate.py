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
    fa_sequence, TR_sequence = np.zeros(info.N_pe + prep_num), (np.ones(info.N_pe + prep_num) * info.TR)
    fa_sequence[0], TR_sequence[0] = info.fa / 2, info.TR / 2
    for i in range (1, info.N_pe + + prep_num):
        fa_sequence[i] = sign((i+1) % 2) * info.fa
    return fa_sequence, TR_sequence
    
def TR_node(fa_sequence, TR_sequence, i, info):
    fa = fa_sequence.tolist()[i]
    t = TR_sequence.tolist()[i]
    init_pulse = {
        "type": "PULSE",
        "FA": fa,
        "t": 0,
        "slice_thickness": info.thickness # 选层
    }
    res = [init_pulse]
    if i == 0: # 只给了一个 prep
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
    
def yaml_generate(info:main_info.info):
    # bssfp
    prep_num = 1
    fa_sequence, TR_sequence = get_sequence_info(info, prep_num)

    seq_list = []
    for i in range(info.N_pe + 1):
        seq_list += TR_node(fa_sequence, TR_sequence, i, info)
    print(seq_list[-5:-1])
    with open("test.yaml", "w") as f:
        yaml.safe_dump(seq_list, f)

yaml_generate(test_info)
# print(sign((0 - 1) % 2) * ((0 + 1) // 2))