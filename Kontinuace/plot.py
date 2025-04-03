import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("C:\Fortran projects\MANDM\Kontinuace\output.txt")

X = data[:, 0]
theta = data[:, 1]
Da = data[:, 2]

plt.figure(figsize=(10, 6))
plt.plot(Da, X, lw='2')

plt.xlabel('Damckeler')
plt.ylabel('Conversion')
plt.grid(True)

plt.show()