import numpy as np
from itertools import product, combinations
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# set the number of monomers and the bond length
N = 100
l = 1.0

# generate random angles
theta = np.random.uniform(0, np.pi, N-2)
phi = np.random.uniform(0, 2*np.pi, N-2)

# compute the positions of the monomers
x = np.zeros(N)
y = np.zeros(N)
z = np.zeros(N)

for i in range(1, N):
    if i == 1:
        x[i] = l * np.sin(theta[i-1]) * np.cos(phi[i-1])
        y[i] = l * np.sin(theta[i-1]) * np.sin(phi[i-1])
        z[i] = l * np.cos(theta[i-1])
    else:
        x[i] = x[i-1] + l * np.sin(theta[i-2]) * np.cos(phi[i-2])
        y[i] = y[i-1] + l * np.sin(theta[i-2]) * np.sin(phi[i-2])
        z[i] = z[i-1] + l * np.cos(theta[i-2])

# center the polymer coordinates
x = x - np.mean(x)
y = y - np.mean(y)
z = z - np.mean(z)

# plot the polymer in 3D with the ball and stick model
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, 'o-', markersize=0.5, linewidth=0.1, color='blue')
for i in range(N):
    ax.scatter(x[i], y[i], z[i], s=40, c='blue', alpha=0.5)
    if i < N-1:
        ax.plot([x[i], x[i+1]], [y[i], y[i+1]], [z[i], z[i+1]], linewidth=1.5, color='black')

# plot the cube with edges
r = [-5, 5]
for s, e in combinations(np.array(list(product(r,r,r))), 2):
    if np.sum(np.abs(s-e)) == r[1]-r[0]:
        ax.plot3D(*zip(s, e), color="black", linewidth=0.5)

# turn off the grid
ax.grid(False)

# turn off the axis labels, ticks, and tick labels
ax.set_axis_off()
plt.tight_layout()
plt.show()

