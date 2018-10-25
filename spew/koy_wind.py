import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

# rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)


def cloud(x, y, v, t, h):
    V = (np.pi * ((x - v * t)**2 + y**2) * h) / t
    return V


v = 30
h = 5000
Vdot = 5 * (10**10)
t = np.linspace(1, 6000, 10)


y = np.linspace(-150000, 150000, 400)
x = np.linspace(-50000, 350000, 400)
x, y = np.meshgrid(x, y)
fig, ax = plt.subplots()
ax.axis('equal')
ax.plot(0, 0, 'r^', ms=6)

for ti in t:
    ax.contour(x, y, cloud(x, y, v, ti, h), [Vdot], colors='k')
plt.title(r'Umbrella Cloud Migration under Regional Wind Effects', fontsize=16)
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$y$ (m)')
plt.tight_layout()
plt.savefig("wind_profile.eps", dpi=200, format='eps')
