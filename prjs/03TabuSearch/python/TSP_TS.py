import random

def distance(node1, node2):
    """
    几何直线距离
    """
    return ((node1[0] - node2[0])**2 + (node1[1] - node2[1])**2)**0.5

def total_distance(route, nodes):
    """
    计算路径的总距离
    传参:nodes = [(x1, y1), (x2, y2), ..., (xn, yn)]
    """
    total = 0
    for i in range(len(route)):
        total += distance(nodes[route[i]], nodes[route[(i+1) % len(route)]])
    return total

def tabu_search(nodes, max_iter=1000, tabu_size=10, aspiration_level=0, frequency_table=None):
    """
    禁忌搜索算法
    参数:
        nodes (list): 包含各个节点坐标的列表
        max_iter (int): 最大迭代次数,默认为1000
        tabu_size (int): 禁忌表的大小,默认为10
        aspiration_level (float): 渴望水平,默认为0,当发现更好的解且优于该水平时,允许违反禁忌策略
        frequency_table (dict): 频数表,记录每个动作的执行次数,用于调整禁忌表的移除和选择
        penalty (float): 惩罚因子取值一般应远小于目标值(1%目标值或1‰目标值),α越大分散性越好,广域搜索能力强,但会损坏邻域搜索

    返回:
        tuple: (best_route, min_distance),最佳路径和最佳路径的总距离
    """
    current_route = list(range(len(nodes)))
    best_route = current_route[:]
    tabu_list = []
    iter_count = 0
    penalty = 0.1
    frequency_table = {}# 初始化频数表

    while iter_count < max_iter:
        best_neighbor = None
        best_distance = float('inf')# 使得范围小于inf

        # 遍历计算每个邻域的数值
        for i in range(len(current_route)):
            for j in range(i + 1, len(current_route)):
                new_route = current_route[:]
                new_route[i], new_route[j] = new_route[j], new_route[i]
                new_distance = total_distance(new_route, nodes)

                # 计算当前动作的频数
                current_frequency = frequency_table.get((i, j), 0) # 获取键,默认0

                # 算后的筛选并不能剪枝来减少运行次数
                if new_distance < best_distance and (i, j) not in tabu_list:
                    # 考虑频数,对动作进行加权
                    weighted_distance = new_distance + current_frequency * penalty
                    if weighted_distance < best_distance:
                        best_neighbor = new_route[:]
                        best_distance = weighted_distance
                        best_i, best_j = i, j


        if best_neighbor is None:
            break

        current_route = best_neighbor[:]
        tabu_list.append((best_i, best_j))
        if len(tabu_list) > tabu_size:
            tabu_list.pop(0)

        if best_distance < total_distance(best_route, nodes):
            best_route = current_route[:]

        # 渴望水平判断
        if best_distance < aspiration_level:
            tabu_list.clear()  # 清空禁忌表

        # 更新频数表
        if frequency_table is not None:
            frequency_table[(best_i, best_j)] = frequency_table.get((best_i, best_j), 0) + 1

        iter_count += 1

    return best_route, total_distance(best_route, nodes)

# 测试
cities = [(1304, 2312), (3639, 1315), (4177, 2244), (3712, 1399), (3488, 1535),
          (3326, 1556), (3238, 1229), (4196, 1044), (4312, 790), (4386, 570),
          (3007, 1970), (2562, 1756), (2788, 1491), (2381, 1676), (1332, 695),
          (3715, 1678), (3918, 2179), (4061, 2370), (3780, 2212), (3676, 2578),
          (4029, 2838), (4263, 2931), (3429, 1908), (3507, 2376), (3394, 2643),
          (3439, 3201), (2935, 3240), (3140, 3550), (2545, 2357), (2778, 2826),
          (2370, 2975)]



best_route, min_distance = tabu_search(cities, max_iter=1000, tabu_size=10, aspiration_level=500)
print("Best Route:", best_route)
print("Minimum Distance:", min_distance)
