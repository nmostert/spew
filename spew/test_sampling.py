import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from eruption import Eruption, plot_grid


filename = './data/ceroNegro_GITHUB.txt'

cn = Eruption(data=filename, vent=Point(532290, 1382690), test=False)
cn
filters = {
    'radius': (5000, np.inf)}
sample = cn.sample(50, weight='MassArea', alpha=0.1,
                   beta=0.1, filters=filters)
fig, ax = plot_grid(sample, vent=cn.vent)
fig, ax = plot_grid(cn.df, vent=cn.vent)
temp_df = cn.df.copy()
for col in filters:
    temp_df = temp_df[(temp_df[col] > filters[col][0]) &
                      (temp_df[col] < filters[col][1])]
norm_mass = temp_df["MassArea"].values / \
    np.max(temp_df["MassArea"].values)


def f1(w, a, b):
    return w ** a


def f2(w, a, b):
    return np.log(w + 1) ** (b * (a ** 2))


def f3(w, a, b):
    return w ** (b * np.log(1 + a))


def f4(w, a, b):
    return w ** (b * (a ** 2))


xx = np.linspace(0, 1, num=200)

a1s = np.linspace(0, 1, num=5)
a2s = np.linspace(0, 1, num=5)
a3s = np.linspace(0, 1, num=5)
a4s = np.linspace(0, 1, num=5)

yf1 = [f1(xx, a, 0) for a in a1s]
yf2 = [f2(xx, a, 0) for a in a2s]
yf3 = [f3(xx, a, 0.1) for a in a3s]
yf4 = [f4(xx, a, 0.1) for a in a4s]

for y1, a1 in zip(yf1, a1s):
    plt.plot(xx, y1)
    plt.legend([r'$\alpha = ' + str(a) + r'$' for a in a1s])
    plt.title(r'$f_1(w) = w^{\alpha}$', fontsize=16)
    plt.xlabel(r'$w$ (Normalised)')
    plt.ylabel(r'$f_1(w)$')
    plt.savefig('f1.eps', dpi=200, format='eps')

for y2, a2 in zip(yf2, a2s):
    plt.plot(xx, y2)
    plt.legend([r'$\alpha = ' + str(a) + r'$' for a in a2s])
    plt.title(r'$f_2(w) = e^{\alpha \sqrt{w}}$', fontsize=16)
    plt.xlabel(r'$w$ (Normalised)')
    plt.ylabel(r'$f_2(w)$')
    plt.savefig('f2.eps', dpi=200, format='eps')

for y3, a3 in zip(yf3, a3s):
    plt.plot(xx, y3)
    plt.legend([r'$\alpha = ' + str(a) + r'$' for a in a2s])
    plt.title(r'$f_3(w) = w^{\beta \log(1+\alpha)}, \quad \beta=0.1$',
              fontsize=16)
    plt.xlabel(r'$w$ (Normalised)')
    plt.ylabel(r'$f_3(w)$')
    plt.savefig('f3.eps', dpi=200, format='eps')

for y4, a4 in zip(yf4, a4s):
    plt.plot(xx, y4)
    plt.legend([r'$\alpha = ' + str(a) + r'$' for a in a4s])
    plt.title(
        r'$f_4(w) = w^{\beta \alpha^2}, \quad \beta=0.1$', fontsize=16)
    plt.xlabel(r'$w$ (Normalised)')
    plt.ylabel(r'$f_4(w)$')
    plt.savefig('f4.eps', dpi=200, format='eps')

tits = [
    r'$f_1(w) = w^{\alpha}$',
    r'$f_2(w) = e^{\alpha \sqrt{w}}$',
    r'$f_3(w) = w^{\beta \log(1+\alpha)}, \quad \beta=0.1$',
    r'$f_4(w) = w^{\beta \alpha^2}, \quad \beta=0.1$'
]
ylabs = [
    r'$f_1(w)$',
    r'$f_2(w)$',
    r'$f_3(w)$',
    r'$f_4(w)$'
]
count = 0
count2 = 0
for f, tit, yl in zip([f1, f2, f3, f4], tits, ylabs):
    count += 1
    count2 = 0
    for a, b in zip([0, 0.1, 0.5, 0.9, 1], [0.1] * 5):
        count2 += 1
        y = f(norm_mass, a, b)
        yr = [yy * np.random.random_sample() for yy in y]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3))
        sc = ax1.scatter(temp_df["Easting"], temp_df["Northing"],
                         marker=',', c=y, vmin=0, vmax=1)
        cbar = fig.colorbar(sc, ax=ax1)

        ax1.plot(cn.vent.x, cn.vent.y, 'r^', ms=8)
        samp = cn.sample(50, weight='MassArea', alpha=a,
                         beta=b, filters=filters, scorefun=f)
        ax1.scatter(samp["Easting"], samp["Northing"],
                    marker='.', color='r', s=5)
        ax2.hist(norm_mass, bins=100, log=True)
        ax3 = ax2.twinx()
        ax3.plot(xx, f(xx, a, b), 'r-')
        ax1.set(aspect='equal')
        ax1.set_title(tit, fontsize=14)
        ax1.set_ylabel(r'Northing')
        ax1.set_xlabel(r'Easting')
        ax2.set_title(r'$\alpha = $ ' + str(a), fontsize=14)
        ax2.set_xlabel(r'$w$ (Normalised)')
        ax2.set_ylabel(r'log(Frequency)')
        ax3.set_ylabel(yl)
        ax3.set_ylim([0, 1.1])
        plt.tight_layout()
        # plt.savefig('samp_f' + str(count) + '_a' + str(count2) +
        # '.eps', dpi=200, format='eps')
