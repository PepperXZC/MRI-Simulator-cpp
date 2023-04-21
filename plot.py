import torch
import numpy as np
import matplotlib.pyplot as plt
# with open("real_mat.txt") as f:
def read(path):
    f = open(path)
    lines = f.readlines()
    data = torch.zeros((64,64), dtype=float)
    for l in range(len(lines)):
        le = lines[l].strip('\n').split(' ')
        r = 0
        for t in le:
            if t != '':
                data[l, r] = float(t)
                r += 1
    return data
real, imag = read("real_mat.txt"), read("img_mat.txt")
data = torch.complex(real, imag)
plt.imshow(data.abs().numpy())
plt.show()
