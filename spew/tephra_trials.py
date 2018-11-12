from eruption import Eruption
from scipy.optimize import curve_fit
import numpy as np
from shapely.geometry import Point
import matplotlib.pyplot as plt


vals = [0, 288, 3600, 7200, 10800, 14400, 18000, 10000000]

param = "fall_time_threshold"

trial = 2

eruptions = dict()

for v in vals:
    filename = './data/%s_trial_%d/%s_trial_%d.txt' % (
        param, trial, param, v)
    eruptions[v] = Eruption(
        data=filename, vent=Point(532290, 1382690), test=False)
    eruptions[v].df = eruptions[v].df[(eruptions[v].df["MassArea"] != 0)
                                      & (eruptions[v].df["MassArea"] > 0.1)
                                      & (eruptions[v].df["MassArea"] < 40000)]
phis = eruptions[vals[0]].phi_labels

rows = int(len(phis) / 2)

for h, erp in eruptions.items():
    xx = erp.df.radius
    yy = erp.df.MassArea
    plt.semilogy(xx, yy, label=r'FTT = %d s' % h, lw=1)
    plt.legend()
    plt.title(
        r'Mass/Area against distance from the vent')
    plt.ylabel('Mass/Area')
    plt.xlabel(r'Radius')
    plt.tight_layout()
# plt.savefig("log_sqr.eps", dpi=200, format='eps')


fig, axs = plt.subplots(rows, 2, figsize=(
    10, 14), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, phi in enumerate(phis):
    for h, erp in eruptions.items():
        new_df = erp.df[erp.df[phi] != 0]
        xx = (new_df.radius)**2
        yy = new_df.MassArea * (new_df[phi] / 100)
        axs[i].semilogy(xx, yy, '-', label=r'FTT = %d s' % h, lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'Log(Mass/Area) against $\mathrm{Radius}^2$ for $\phi$ class %s' % phi)
        axs[i].set_ylabel('Log(Mass/Area)')
        axs[i].set_xlabel(r'Radius$^2$')
        # axs[i].set_xlim(0, 1.7e10)
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/koy_plot_trial_%d.pdf" %
            (param, trial, trial), dpi=200, format='pdf')

fig, axs = plt.subplots(rows, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, phi in enumerate(phis):
    for h, erp in eruptions.items():
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df[phi]
        axs[i].semilogx(xx, yy, '-', ms=4, label=r'FTT = %d s' % h, lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'Weight percentage of $\phi$ class %s against $\mathrm{Radius}$' % phi)
        axs[i].set_ylabel(r'Weight percentage for $\phi \in %s$' % phi)
        axs[i].set_xlabel(r'Radius')
        # axs[i].set_ylim(top=100)
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/phi_over_rad_plot_trial_%d.pdf" %
            (param, trial, trial), dpi=200, format='pdf')

fig, axs = plt.subplots(4, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, (h, erp) in enumerate(eruptions.items()):
    for phi in phis:
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df[phi]
        axs[i].semilogx(xx, yy, '-', ms=4, label=r'$\phi = %s$' % phi, lw=1)
        axs[i].legend()  # bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        axs[i].set_title(
            r'Weight percentage against $\mathrm{Radius}$, %d s FTT' % h)
        axs[i].set_ylabel(r'Weight percentage')
        axs[i].set_xlabel(r'Radius')
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/param_over_rad_plot_trial_%d.pdf" %
            (param, trial, trial), dpi=200, format='pdf')


def func(x, a, b, c, d):
    return a * (b ** (c * x)) + d


def func2(x, a, b, c):
    return a * np.exp(-b * x) + c


def func3(x, m, c):
    return m * x + c


fig, axs = plt.subplots(3, 2, figsize=(
    10, 12), facecolor='w', edgecolor='k')
axs = axs.ravel()
slopes = []
intercepts = []
for i, phi in enumerate(phis):
    m = []
    c = []
    for h, erp in eruptions.items():

        new_df = erp.df[erp.df[phi] != 0]
        xx = new_df.radius**2
        yy = np.log(new_df.MassArea * (new_df[phi] / 100))
        axs[i].plot(xx, yy, 'k*', ms=3, label=None, alpha=0.5)

        popt, pcov = curve_fit(func3, xx, yy)
        m.append(popt[0])
        c.append(popt[1])
        axs[i].plot(xx, func3(xx, *popt), '-',
                    label='H = %d, m=%.2e, c=%5.2f' % (h, popt[0], popt[1]))
    axs[i].legend()
    axs[i].set_ylabel('Log(Mass/Area)')
    axs[i].set_xlabel(r'Radius$^2$')
    axs[i].set_title(r'$\phi \in %s$' %
                     phi, fontsize=14)
    # axs[i].set_xlim(0, 1.7e10)
    slopes.append(m)
    intercepts.append(c)
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/line_fits_trial_%d.png" %
            (param, trial, trial), dpi=200, format='png')

for phi, slope in zip(phis, slopes):
    plt.plot(vals, slope, label=phi, marker='o')
plt.title(r'Plot of slope against column height', fontsize=14)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.xlabel("Column Height (km)")
plt.ylabel("Slope")
plt.tight_layout()
plt.legend()
# plt.show()
plt.savefig("./data/%s_trial_%d/slope_trial_%d.png" %
            (param, trial, trial), dpi=200, format='png')

for phi, intercept in zip(phis, intercepts):
    plt.plot(vals, intercept, label=r'\(\phi \in %s\)' %
             phi, marker='o')
plt.title(r'Plot of intercept against column height', fontsize=14)
plt.xlabel("Column Height (km)")
plt.ylabel("Intercept")
plt.tight_layout()
plt.legend()
# plt.show()
plt.savefig("./data/%s_trial_%d/intercept_trial_%d.png" %
            (param, trial, trial), dpi=200, format='png')
