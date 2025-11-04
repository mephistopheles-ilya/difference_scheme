import numpy as np
import matplotlib.pyplot as plt

file_name = "2d.txt"
dx = 0.001
dt = 0.01



dat = []
with open(file_name, "r") as file:
    for i, line in enumerate(file):
        elements = line.split()
        row = [float(elements[i]) for i in range(0, len(elements), 1)]
        dat.append(row)

H_data = np.array(dat)
dat.clear()

plt.figure(figsize=(10, 8))
plt.imshow(H_data.T, aspect='auto', cmap='hot', extent=(0, 216.3, 0, 10))

plt.colorbar(label='H(x, t)')
plt.xlabel('Time')
plt.ylabel('Space')
plt.show()

'''
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
plt.ylabel('H(x)')
plt.legend()
plt.show()
'''
