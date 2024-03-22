import numpy as np
import matplotlib.pyplot as plt

# 提供的数据
x = np.array([150, 125, 75, 50])
y = np.array([125, 122, 113, 104])

# 拟合曲线
p = np.polyfit(np.log(x), y, 2)
x_fit = np.linspace(0.1, 4, 1000)
y_fit = np.polyval(p, np.log(x_fit))

# 绘制曲线
plt.plot(x_fit, y_fit)

# 添加标题和标签
plt.title('Smooth Curve Simulation for profile time')
plt.xlabel('control variable')
plt.ylabel('profile time')

# 显示图形
plt.show()
