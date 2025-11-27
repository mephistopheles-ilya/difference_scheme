import numpy as np
import matplotlib.pyplot as plt

def H(t, x):
    return np.exp(t) * (np.cos(3 * np.pi * x) + 1.5)

def V(t, x):
    return np.cos(2 * np.pi * t) * np.sin(4 * np.pi * x)

h_a = 0.0001
h_c = 0.01
L = 1.0
x_a = np.arange(0, L + h_a, h_a) 
x_c = np.arange(0, L + h_c, h_c) 

# Create a figure with 2 subplots side by side
plt.figure(figsize=(12, 5))

# First subplot - H
plt.subplot(1, 2, 1)
filename = "H_res.txt"
H_c = np.loadtxt(filename)
H_a = H(1, x_a)

plt.plot(x_c, H_c, label="calc")
plt.plot(x_a, H_a, label="real")
plt.legend()
plt.title("H function")
plt.xlabel("x")
plt.ylabel("H(1, x)")

# Second subplot - V
plt.subplot(1, 2, 2)
filename = "V_res.txt"
V_c = np.loadtxt(filename)
V_a = V(1, x_a)

plt.plot(x_c, V_c, label="calc")
plt.plot(x_a, V_a, label="real")
plt.legend()
plt.title("V function")
plt.xlabel("x")
plt.ylabel("V(1 , x)")

plt.tight_layout()
plt.show()
