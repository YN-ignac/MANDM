import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("C:\Fortran projects\MANDM\Bifurcation-diagram\output.txt")

X = data[:, 0]
theta = data[:, 1]
Da = data[:, 2]
beta = data[:, 3]

plt.figure(figsize=(10, 6))
plt.plot(Da, beta, lw='2')

plt.xlabel('Damckeler')
plt.ylabel('Beta')
plt.grid(True)

plt.show()