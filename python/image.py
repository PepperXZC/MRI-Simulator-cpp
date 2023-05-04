def get_point_index(length, bandwidth):
    # 矩形
    # TODO：生成两个 point_index 集合
    li_vassel, li_muscle = [], []
    a, b = int(length // 2), int(bandwidth // 2)
    lower, upper = a - b, a + b
    for i in range(length):
        for j in range(length):
            if j >= lower + 1 and j < upper + 1:
                li_vassel.append((i, j))
            else:
                li_muscle.append((i, j))
    return li_vassel, li_muscle
