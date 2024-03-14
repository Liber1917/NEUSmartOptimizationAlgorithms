import numpy as np
import random
import matplotlib.pyplot as plt

# 计算路径总长度（根据城市坐标）
def func1(C, solution):
    total_distance = 0
    for i in range(len(solution) - 1):
        city1 = solution[i]
        city2 = solution[i + 1]
        total_distance += np.sqrt((C[city1, 0] - C[city2, 0]) ** 2 + (C[city1, 1] - C[city2, 1]) ** 2)
    # 添加回到起始城市的距离
    total_distance += np.sqrt((C[solution[-1], 0] - C[solution[0], 0]) ** 2 + (C[solution[-1], 1] - C[solution[0], 1]) ** 2)
    return total_distance

# 城市坐标
C = np.array([[1304, 2312], [3639, 1315], [4177, 2244], [3712, 1399], [3488, 1535], [3326, 1556],
              [3238, 1229], [4196, 1044], [4312, 790], [4386, 570], [3007, 1970], [2562, 1756],
              [2788, 1491], [2381, 1676], [1332, 695], [3715, 1678], [3918, 2179], [4061, 2370],
              [3780, 2212], [3676, 2578], [4029, 2838], [4263, 2931], [3429, 1908], [3507, 2376],
              [3394, 2643], [3439, 3201], [2935, 3240], [3140, 3550], [2545, 2357], [2778, 2826],
              [2370, 2975]])
N = len(C)  # 城市数目

Tabu = np.zeros((N, N))  # 禁忌表
TabuL = int(np.sqrt((N * (N - 1)) / 2))  # 禁忌长度
Ca = 200  # 领域解个数
CaNum = np.zeros((Ca, N), dtype=int)  # 领域解种群
S0 = np.random.permutation(N)  # 随机产生初始解
bestsofar = S0.copy()  # 当前最佳解
BestL = np.inf  # 当前最佳解距离
Gmax = 1000  # 最大迭代次数

# 初始化禁忌搜索循环
p = 0
ALong = np.zeros(Gmax)
ArrBestL = np.zeros(Gmax)

# 禁忌搜索循环
while p < Gmax:
    ALong[p] = func1(C, S0)  # 当前解适配值
    i = 0
    A = np.zeros((Ca, 2), dtype=int)  # 交换的城市矩阵
    while i < Ca:
        M = np.random.randint(N, size=2)
        if M[0] != M[1]:
            A[i, 0] = max(M)
            A[i, 1] = min(M)
            if i == 0:
                isa = 0
            else:
                isa = np.any((A[i, 0] == A[:i, 0]) & (A[i, 1] == A[:i, 1]))
            if not isa:
                i += 1
        else:
            pass
    F = np.zeros(Ca)
    for i in range(Ca):
        CaNum[i] = S0
        CaNum[i, [A[i, 1], A[i, 0]]] = S0[[A[i, 0], A[i, 1]]]
        F[i] = func1(C, CaNum[i])
    BestCaNum = Ca // 2
    BestCa = np.full((BestCaNum, 4), np.inf)
    for i in range(Ca):
        if i < BestCaNum:
            BestCa[i, 0] = i
            BestCa[i, 1] = F[i]
            BestCa[i, 2] = S0[A[i, 0]]
            BestCa[i, 3] = S0[A[i, 1]]
        else:
            idx = np.argmin(BestCa[:, 1])
            if F[i] < BestCa[idx, 1]:
                BestCa[idx, 0] = i
                BestCa[idx, 1] = F[i]
                BestCa[idx, 2] = S0[A[i, 0]]
                BestCa[idx, 3] = S0[A[i, 1]]
    BestCa = BestCa[np.argsort(BestCa[:, 1])]
    if BestCa[0, 1] < BestL:
        BestL = BestCa[0, 1]
        S0 = CaNum[int(BestCa[0, 0])]
        bestsofar = S0.copy()
        Tabu -= 1
        Tabu[int(A[0, 1]), int(A[0, 0])] = TabuL
    else:
        for i in range(BestCaNum):
            if Tabu[int(BestCa[i, 2]), int(BestCa[i, 3])] == 0:
                S0 = CaNum[int(BestCa[i, 0])]
                Tabu -= 1
                Tabu[int(BestCa[i, 2]), int(BestCa[i, 3])] = TabuL
                break
    ArrBestL[p] = BestL
    p += 1

# 绘制适应度进化曲线
plt.figure(2)
plt.plot(ArrBestL)
plt.xlabel('迭代次数')
plt.ylabel('目标函数值')
plt.title('适应度进化曲线')
plt.show()

# 输出结果
BestShortcut = bestsofar  # 最佳路线
theMinDistance = BestL  # 最佳路线长度
print("最佳路径：", BestShortcut)
print("最短距离：", theMinDistance)
