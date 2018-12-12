from eruption import Eruption
from scipy.optimize import curve_fit, minimize
import numpy as np
from shapely.geometry import Point
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import scipy.ndimage
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc


rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)

mpl.style.use("ggplot")


param = "plume_height"
label = "H"
units = "km"
trial = 7
vals = np.linspace(10000, 40000, 21)

vals = [int(val) for val in vals]


def disp_func(col):
    return col / 1000


eruptions = dict()

for v in vals:
    filename = './data/%s_trial_%d/%s_trial_%d.txt' % (
        param, trial, param, v)
    eruptions[v] = Eruption(
        data=filename, vent=Point(532290, 1382690), test=False)
    eruptions[v].df = eruptions[v].df[(eruptions[v].df["MassArea"] != 0)
                                      & (eruptions[v].df["MassArea"] > 0.1)
                                      & (eruptions[v].df["radius"] < 40000)]

for key in eruptions:
    print(key)

phis = eruptions[vals[0]].phi_labels

display(eruptions[40000].df)
rows = int(len(phis) / 2)
rows


def gamma_kernel(x, a, b, c):
    return c * (x**a) * np.exp(b * x)


def sse(k, x_data, y_data):
    return sum((y_data - (gamma_kernel(x_data, *k))) ** 2)


fit_df = pd.DataFrame(columns=['H', 'Phi', 'a', 'b', 'c'])

for k, phi in enumerate(phis):
    # fig, axs = plt.subplots(rows, 2, figsize=(
    #     10, 15), facecolor='w', edgecolor='k')
    # axs = axs.ravel()
    for i, (h, erp) in enumerate(eruptions.items()):
        new_df = erp.df
        ma = new_df['MassArea'].values
        wt = new_df[phi].values
        xx = new_df['radius'].values
        yy = ma * (wt / 100)

        pnts = list(set(zip(xx, yy)))
        pnts = sorted(pnts, key=lambda x: x[0])

        x = np.array([p[0] for p in pnts], dtype=np.float64)
        y = np.array([p[1] for p in pnts], dtype=np.float64)

        tot_mass = np.trapz(y, x=x)
        idx = 0
        perc = 0.0
        while perc < 99.99999999:
            mass = np.trapz(y[:idx], x=x[:idx])
            idx += 1
            perc = (mass / tot_mass) * 100

        x_data = x[:idx]
        y_data = y[:idx]

        type(x_data)

        k0 = np.array([0.0, 0.0, 0.0], dtype=np.float64)

        def fun(k): return sse(k, x_data, y_data)

        popt = minimize(fun, k0)
        display(popt)
        xf = np.linspace(x_data.min(), x_data.max(), 100)
        yf = gamma_kernel(xf, *popt.x)

        plt.scatter(x_data, y_data, s=30, facecolors='none', edgecolors='r')
        plt.plot(xf, yf, '--', color='k', lw=1.2, label=r'\noindent $a$=%.2e,\\ $b$=%.2e,\\ $c$=%.2e' % (
            popt.x[0], popt.x[1], popt.x[2]))
        # axs[i].legend()
        # axs[i].set_ylabel('Mass/Area')
        # axs[i].set_xlabel(r'Radius')
        plt.title(r'$\phi \in %s$, $H = %d$km' %
                  (phi, disp_func(h)), fontsize=14)
        fit_df = fit_df.append(
            {'H': h, 'Phi': phi, 'a': popt.x[0], 'b': popt.x[1], 'c': popt.x[2]}, ignore_index=True)
        plt.show()
    # plt.tight_layout()
    # plt.show()
    # plt.savefig("./data/%s_trial_%d/minimize_%d.pdf" %
    #             (param, trial, k), format='pdf')

fit_df

a_piv = fit_df.pivot_table(index='Phi', columns='H', values='a')
b_piv = fit_df.pivot_table(index='Phi', columns='H', values='b')
c_piv = fit_df.pivot_table(index='Phi', columns='H', values='c')
a_piv = a_piv.reindex(phis)
b_piv = b_piv.reindex(phis)
c_piv = c_piv.reindex(phis)
phis
a_piv

xlabs = a_piv.columns.values / 1000
xlabs
xlabels = [123] + list(xlabs)
ylabels = ['213'] + phis

fig, axs = plt.subplots(3, 1, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, piv, p in zip([0, 1, 2], [a_piv, b_piv, c_piv], ['a', 'b', 'c']):
    c = axs[i].imshow(piv.values, cmap='inferno')
    # cont = axs[i].contour(scipy.ndimage.zoom(piv.values, 3), colors='w')
    axs[i].set_xticklabels(xlabels)
    axs[i].set_yticklabels(ylabels)
    fig.colorbar(c, ax=axs[i])
    axs[i].grid(False)
    axs[i].set_title(p)
    axs[i].clabel(cont, inline=1, fontsize=10)
    axs[i].set_xlabel("Column height (km)")
    axs[i].set_ylabel("Phi class")
# plt.show()
plt.savefig("./data/%s_trial_%d/param_maps_cont.pdf" %
            (param, trial), format='pdf')
