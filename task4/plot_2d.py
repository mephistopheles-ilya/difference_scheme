import numpy as np
import matplotlib.pyplot as plt

file_name = "2d.txt"
dx = 0.001



dat = []
with open(file_name, "r") as file:
    for i, line in enumerate(file):
        if (i % 2 == 0):
            elements = line.split()
            row = [float(elements[i]) for i in range(0, len(elements), 1)]
            dat.append(row)

dat.pop()
H_data = np.array(dat)
dat.clear()

stab_index = H_data.shape[0]
y_values1 = H_data[int(stab_index / 4)]
y_values2 = H_data[int(stab_index / 2)]
y_values3 = H_data[int(2 * stab_index / 3)]
y_values4 = H_data[int(stab_index - 1)]
x_values = np.arange(len(y_values1)) * dx

plt.figure(figsize=(10, 6))
plt.plot(x_values, y_values1, label='n / 4')
plt.plot(x_values, y_values2, label='n / 2')
plt.plot(x_values, y_values3, label='2 * n / 3')
plt.plot(x_values, y_values4, label='n')
plt.xlabel('x')
plt.ylabel('V(x)')
plt.legend()
plt.show()


