import numpy as np

# 种群数量
pop_size = 10
# 交叉概率
PC = 0.6
# 变异概率
PM = 0.01
# 最大值
X_max = 5
# 最小值
X_min = 0
# DNA长度与保留位数有关, 2^10 当前保留3位小数点
DNA_SIZE = 10
# 进化的代数
N_GENERATIONS = 1000

# 目标函数，根据给定的x值计算y值
def aim(x):
    return 10 * np.sin(5 * x) + 7 * np.cos(4 * x)

# 解码函数，将二进制编码转换为实数值
def decode(pop):
    return pop.dot(2 ** np.arange(DNA_SIZE)[::-1]) * (X_max - X_min) / float(2**DNA_SIZE - 1) + X_min

# 适应度函数
def fitnessget(pred):
    return pred + 1e-3 - np.min(pred) # 相对差异

# 选择函数，根据适应度进行选择
def select(pop, fitness):
    idx = np.random.choice(np.arange(pop_size), size=pop_size, replace=True, p=fitness/fitness.sum()) # 根据适应度进行随机选择，适应度高的被选中的概率大
    return pop[idx]

# 交叉操作
def change(parent, pop):
    if np.random.rand() < PC: # 随机数小于某个交叉概率PC，则执行交叉操作
        i_ = np.random.randint(0, pop_size, size=1) # 随机选择另一个个体作为交叉对象
        cross_points = np.random.randint(0, 2, size=DNA_SIZE).astype(bool) # 随机生成数组，元素表示交叉点位置
        parent[cross_points] = pop[i_, cross_points] # 基选中因复制到当前个体，实现交叉
    return parent

# 变异操作
def variation(child, pm):
    for point in range(DNA_SIZE):
        if np.random.rand() < pm:
            child[point] = 1 if child[point] == 0 else 0 # 随机数小于某个变异概率pm，则执行变异操作(取反)
    return child

# 初始化种群
pop = np.random.randint(2, size=(pop_size, DNA_SIZE))

# 开始迭代进化
for i in range(N_GENERATIONS):
    # 解码
    X_value = decode(pop)

    # 获取目标函数值
    F_values = aim(X_value)

    # 获取适应值
    fitness = fitnessget(F_values)

    # 记录最大值及对应的DNA
    if i == 0: # 第一代直接记录
        max = np.max(F_values)
        max_DNA = pop[np.argmax(F_values), :]

    if max < np.max(F_values):
        max = np.max(F_values)
        max_DNA = pop[np.argmax(F_values), :]

    # 每10代打印最适应值和对应的X值
    if i % 10 == 0:
        print("Most fitted value and X: \n", np.max(F_values), decode(pop[np.argmax(F_values), :]))

    # 选择
    pop = select(pop, fitness)

    # 复制种群
    pop_copy = pop.copy()

    # 对每个父代进行交叉和变异
    for parent in pop:
        child = change(parent, pop_copy)
        child = variation(child, PM)
        parent[:] = child

# 输出最大的目标函数值及对应的DNA和X值
print(" maximum value of the objective function :", max)
print("DNA value", max_DNA)
print("X value :", decode(max_DNA))
