from eruption import Eruption
from scipy.optimize import curve_fit
import numpy as np
from shapely.geometry import Point
import matplotlib.pyplot as plt


col_heights = [5, 10, 15, 20, 25, 30, 35, 40]

eruptions = dict()

trial = 5

for h in col_heights:
    if trial == 1:
        filename = './data/col_trial/col_trail_%d.txt' % h
    else:
        filename = './data/col_trial_%d/col_trail_%d.txt' % (trial, h)
    eruptions[h] = Eruption(
        data=filename, vent=Point(532290, 1382690), test=False)
    eruptions[h].df = eruptions[h].df[(eruptions[h].df["MassArea"] != 0)
                                      & (eruptions[h].df["MassArea"] > 0.1)
                                      & (eruptions[h].df["MassArea"] < 40000)]

len(eruptions[40].df)

phis = eruptions[25].phi_labels
#
# for phi in phis:
#     for h, erp in eruptions.items():
#         xx = erp.df.radius
#         yy = erp.df.MassArea * (erp.df[phi] / 100)
#         plt.plot(xx, yy, '.', label=r'H = %d km' % h, lw=1)
#         plt.legend()
#         plt.title(r'Mass/Area against Radius from the vent')
#         plt.ylabel('Mass/Area')
#         plt.xlabel('Radius (m)')
#         plt.tight_layout()
#     plt.show()
# # plt.savefig("lin_lin.eps", dpi=200, format='eps')
#
# for h, erp in eruptions.items():
#     xx = erp.df.radius
#     yy = np.log(erp.df.MassArea)
#     plt.plot(xx, yy, label=r'H = %d km' % h, lw=1)
#     plt.legend()
#     plt.title(r'Log(Mass/Area) against Radius from the vent')
#     plt.ylabel('Log(Mass/Area)')
#     plt.xlabel('Radius (m)')
#     plt.tight_layout()
# # plt.savefig("log_lin.eps", dpi=200, format='eps')
#
# for h, erp in eruptions.items():
#     xx = np.log(erp.df.radius)
#     yy = np.log(erp.df.MassArea)
#     plt.plot(xx, yy, label=r'H = %d km' % h, lw=1)
#     plt.legend()
#     plt.title(r'Log(Mass/Area) against Log(Radius) from the vent')
#     plt.ylabel('Log(Mass/Area)')
#     plt.xlabel('Log(Radius)')
#     plt.tight_layout()
# # plt.savefig("log_log.eps", dpi=200, format='eps')
#
# for h, erp in eruptions.items():
#     xx = np.sqrt(erp.df.radius)
#     yy = np.log(erp.df.MassArea)
#     plt.plot(xx, yy, label=r'H = %d km' % h, lw=1)
#     plt.legend()
#     plt.title(r'Log(Mass/Area) against $\sqrt{\mathrm{Radius}}$ from the vent')
#     plt.ylabel('Log(Mass/Area)')
#     plt.xlabel(r'$\sqrt{\mathrm{Radius}}$')
#     plt.tight_layout()
# # plt.savefig("log_sqrt.eps", dpi=200, format='eps')
#
# for h, erp in eruptions.items():
#     xx = (erp.df.radius)**2
#     yy = np.log(erp.df.MassArea)
#     plt.plot(xx, yy, label=r'H = %d km' % h, lw=1)
#     plt.legend()
#     plt.title(
#         r'Log of normalised Mass/Area against $\mathrm{Radius}^2$ from the vent')
#     plt.ylabel('Log of normalised Mass/Area')
#     plt.xlabel(r'Radius$^2$')
#     plt.tight_layout()
# # plt.savefig("log_sqr.eps", dpi=200, format='eps')


phis = eruptions[25].phi_labels

fig, axs = plt.subplots(4, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, phi in enumerate(phis):
    for h, erp in eruptions.items():
        new_df = erp.df.sort_values("radius")
        xx = (new_df.radius)**2
        yy = np.log(new_df.MassArea * (new_df[phi] / 100))
        axs[i].plot(xx, yy, '.', label=r'H = %d km' % h, lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'Log(Mass/Area) against $\mathrm{Radius}^2$ for $\phi$ class %s' % phi)
        axs[i].set_ylabel('Log(Mass/Area)')
        axs[i].set_xlabel(r'Radius$^2$')
        # axs[i].set_xlim(0, 1.7e10)
plt.tight_layout()
plt.show()
# plt.savefig("./data/col_trial_%d/koy_plot_trial_%d.png" %
#             (trial, trial), dpi=200, format='png')

fig, axs = plt.subplots(4, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, phi in enumerate(phis):
    for h, erp in eruptions.items():
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df[phi]
        axs[i].plot(xx, yy, 'o-', ms=4, label=r'H = %d km' % h, lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'$\phi \in %s$' % phi)
        axs[i].set_ylabel(r'Weight percentage for $\phi \in %s$' % phi)
        axs[i].set_xlabel(r'Radius')
plt.tight_layout()
plt.show()
# plt.savefig("./data/col_trial_%d/phi_over_rad_plot_trial_%d.png" %
#             (trial, trial), dpi=200, format='png')

fig, axs = plt.subplots(4, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, phi in enumerate(phis):
    for h, erp in eruptions.items():
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df.MassArea * (new_df[phi] / 100)
        axs[i].plot(xx, yy, 'o-', ms=1, label=r'H = %d km' % h, lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'$\phi \in %s$' % phi)
        axs[i].set_ylabel(r'Mass/Area for $\phi \in %s$' % phi)
        axs[i].set_xlabel(r'Radius')
        axs[i].set_ylim(top=70000)
plt.tight_layout()
plt.show()
# plt.savefig("./data/col_trial_%d/phi_over_rad_plot_trial_%d.png" %
#             (trial, trial), dpi=200, format='png')

fig, axs = plt.subplots(4, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, (h, erp) in enumerate(eruptions.items()):
    for phi in phis:
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df[phi]
        axs[i].plot(xx, yy, 'o', ms=4, label=r'$\phi = %s$' % phi, lw=1)
        axs[i].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        axs[i].set_title(
            r'Weight percentage against $\mathrm{Radius}$, %d km column height' % h)
        axs[i].set_ylabel(r'Weight percentage')
        axs[i].set_xlabel(r'Radius')
plt.tight_layout()
# plt.show()
plt.savefig("./data/col_trial_%d/col_over_rad_plot_trial_%d.png" %
            (trial, trial), dpi=200, format='png')


fig, axs = plt.subplots(4, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, (h, erp) in enumerate(eruptions.items()):
    for phi in phis:
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df.MassArea * (new_df[phi] / 100)
        axs[i].plot(xx, yy, 'o-', ms=4, label=r'$\phi = %s$' % phi, lw=1)
        axs[i].legend()  # bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        axs[i].set_title(
            r'Weight percentage against $\mathrm{Radius}$, %d km column height' % h)
        axs[i].set_ylabel(r'Weight percentage')
        axs[i].set_xlabel(r'Radius')
plt.tight_layout()
plt.show()
# plt.savefig("./data/col_trial_%d/col_over_rad_plot_trial_%d.png" %
# (trial, trial), dpi=200, format='png')


def func(x, a, b, c, d):
    return a * (b ** (c * x)) + d


def func2(x, a, b, c):
    return a * np.exp(-b * x) + c


def func3(x, m, c):
    return m * x + c


fig, axs = plt.subplots(2, 4, figsize=(
    20, 10), facecolor='w', edgecolor='k')
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
plt.savefig("koy_fits_yo.png", dpi=200, format='png')

for phi, slope in zip(phis, slopes):
    plt.plot(col_heights, slope, label=phi, marker='o')
plt.title(r'Slope', fontsize=14)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.xlabel("Column Height (km)")
plt.ylabel("Slope")
plt.tight_layout()
plt.legend()
# plt.show()
plt.savefig("slope.png", dpi=200, format='png')

for phi, intercept in zip(phis, intercepts):
    plt.plot(col_heights, intercept, label=r'\(\phi \in %s\)' %
             phi, marker='o')
plt.title(r'Intercept', fontsize=14)
plt.xlabel("Column Height (km)")
plt.ylabel("Intercept")
plt.tight_layout()
plt.legend()
# plt.show()
plt.savefig("intercept.png", dpi=200, format='png')
