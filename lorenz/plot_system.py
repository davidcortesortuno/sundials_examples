import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('lorenz.txt')

plt.plot(data[:, 1], data[:, 2])
plt.xlabel('x')
plt.ylabel('y')

plt.savefig('lorenz_attractor.pdf', bbox_inches='tight')
