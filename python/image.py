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


def get_3_point_index(length, bandwidth):
    # 矩形
    # TODO：生成两个 point_index 集合
    center = [15, 55, 105]
    b = int(bandwidth // 2)
    li_vassel, li_muscle = [[] for _ in range(3)], []

    for i in range(length):
        for j in range(length):
            flag = 0
            for c in range(len(center)):
                __center = 128 - center[c]
                lower, upper = __center - b, __center + b
                if j >= lower + 1 and j < upper + 1:
                    li_vassel[c].append((i, j))
                    flag = 1
                    break
            if flag == 0:
                li_muscle.append((i, j))
    return li_vassel, li_muscle


# li_vassel, li_muscle = get_3_point_index(128, 10)
