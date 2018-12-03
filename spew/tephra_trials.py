from eruption import Eruption
from scipy.optimize import curve_fit
import numpy as np
from shapely.geometry import Point
import matplotlib.pyplot as plt


param = "plume_height"
label = "H"
units = "km"
trial = 6

vals = [5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000]


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
                                      & (eruptions[v].df["MassArea"] < 40000)]
    # filename_csv = './data/%s_trial_%d/%s_trial_%d.csv' % (
    #     param, trial, param, v)
    # eruptions[v].df.to_csv(filename_csv, ',')


phis = eruptions[vals[0]].phi_labels

eruptions[40000].df["MassArea"].sum()

rows = int(len(phis) / 2)
rows


for v, erp in eruptions.items():
    xx = erp.df.radius
    yy = erp.df.MassArea
    plt.semilogy(xx, yy, '.-', label=r'%s = %d %s' %
                 (label, disp_func(v), units), lw=1, ms=1)
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
    for v, erp in eruptions.items():
        new_df = erp.df[erp.df[phi] != 0]
        xx = (new_df.radius)**2
        yy = new_df.MassArea * (new_df[phi] / 100)
        axs[i].semilogy(xx, yy, '-', label=r'%s = %d %s' %
                        (label, disp_func(v), units), lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'Log(Mass/Area) against $\mathrm{Radius}^2$ for $\phi$ class %s' % phi)
        axs[i].set_ylabel('Log(Mass/Area)')
        axs[i].set_xlabel(r'Radius$^2$')
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/logmass_rsqr_by_H_sep_phi.png" %
            (param, trial), dpi=200, format='png')

fig, axs = plt.subplots(rows, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, phi in enumerate(phis):
    for v, erp in eruptions.items():
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df[phi]
        axs[i].semilogx(xx, yy, '.-', ms=4, label=r'%s = %d %s' %
                        (label, disp_func(v), units), lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'Weight percentage of $\phi \in %s$ against $\mathrm{Radius}$' % phi)
        axs[i].set_ylabel(r'Weight percentage')
        axs[i].set_xlabel(r'Radius')
        # axs[i].set_ylim(top=100)
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/wt_by_H_sep_phi_logx.png" %
            (param, trial), dpi=200, format='png')

fig, axs = plt.subplots(rows, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, (v, erp) in enumerate(eruptions.items()):
    for phi in phis:
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df[phi]
        axs[i].plot(xx, yy, '.-', ms=4, label=r'$\phi = %s$' % phi, lw=1)
        axs[i].legend()  # bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        axs[i].set_title(
            r'Normalised grain size against distance, %s = %d %s' %
            (label, disp_func(v), units))
        axs[i].set_ylabel(r'Weight percentage')
        axs[i].set_xlabel(r'Radius')
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/wt_by_phi_sep_H.png" %
            (param, trial), dpi=200, format='png')


fig, axs = plt.subplots(4, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, (v, erp) in enumerate(eruptions.items()):
    for phi in phis:
        new_df = erp.df.sort_values("radius")
        xx = new_df.radius
        yy = new_df.MassArea * (new_df[phi] / 100)
        axs[i].semilogy(xx, yy, '.-', ms=4, label=r'$\phi = %s$' % phi, lw=1)
        axs[i].legend()  # bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        axs[i].set_title(
            r'Mass by grain size against distance, %s = %d %s' %
            (label, disp_func(v), units))
        axs[i].set_ylabel(r'Mass/Area')
        axs[i].set_xlabel(r'Radius')
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/mass_by_phi_sep_H_logx.png" %
            (param, trial), dpi=200, format='png')

fig, axs = plt.subplots(rows, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, phi in enumerate(phis):
    for v, erp in eruptions.items():
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df[phi]
        axs[i].semilogx(xx, yy, '.-', ms=4, label=r'%s = %d %s' %
                        (label, disp_func(v), units), lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'$\phi \in %s$' % phi)
        axs[i].set_ylabel(r'Weight percentage')
        axs[i].set_xlabel(r'Radius')
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/wt_by_H_sep_phi_logx.png" %
            (param, trial), dpi=200, format='png')

fig, axs = plt.subplots(rows, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, phi in enumerate(phis):
    for v, erp in eruptions.items():
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius ** 2
        yy = np.log(new_df.MassArea * (new_df[phi] / 100))
        axs[i].plot(xx, yy, '.-', ms=4, label=r'%s = %d %s' %
                    (label, disp_func(v), units), lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'$\phi \in %s$' % phi)
        axs[i].set_ylabel(r'Mass/Area')
        axs[i].set_xlabel(r'Radius')
plt.tight_layout()
plt.show()
# plt.savefig("./data/%s_trial_%d/mass_by_H_sep_phi.png" %
#             (param, trial), dpi=200, format='png')

fig, axs = plt.subplots(rows, 2, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, phi in enumerate(phis):
    for v, erp in eruptions.items():
        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df.MassArea * (new_df[phi] / 100)
        axs[i].semilogx(xx, yy, '.-', ms=4, label=r'%s = %d %s' %
                        (label, disp_func(v), units), lw=1)
        axs[i].legend()
        axs[i].set_title(
            r'$\phi \in %s$' % phi)
        axs[i].set_ylabel(r'Mass/Area')
        axs[i].set_xlabel(r'Radius')
plt.tight_layout()
# plt.show()
plt.savefig("./data/%s_trial_%d/mass_by_H_sep_phi_logx.png" %
            (param, trial), dpi=200, format='png')


def func3(x, m, c):
    return m * x + c


fig, axs = plt.subplots(rows, 2, figsize=(
    10, 12), facecolor='w', edgecolor='k')
axs = axs.ravel()
slopes = []
intercepts = []
for i, phi in enumerate(phis):
    m = []
    c = []
    for h, erp in eruptions.items():

        new_df = erp.df.sort_values("radius")  # [erp.df[phi] != 0]
        xx = new_df.radius
        yy = new_df.MassArea * (new_df[phi] / 100)
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
plt.show()
# plt.savefig("./data/%s_trial_%d/line_fits_trial_%d.png" %
#             (param, trial, trial), dpi=200, format='png')

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


def func4(x, a, b, c, d):
    return a ** (b * (x + c)) + d


fig, axs = plt.subplots(rows, 2, figsize=(
    10, 12), facecolor='w', edgecolor='k')
axs = axs.ravel()
slopes = []
intercepts = []
for i, phi in enumerate(phis):
    a = []
    b = []
    c = []
    d = []
    for h, erp in eruptions.items():
        new_df = erp.df[erp.df[phi] != 0]
        xx = new_df.radius**2
        yy = np.log(new_df.MassArea * (new_df[phi] / 100))
        axs[i].plot(xx, yy, 'k*', ms=3, label=None, alpha=0.5)

        popt, pcov = curve_fit(func3, xx, yy)
        a.append(popt[0])
        b.append(popt[1])
        c.append(popt[2])
        d.append(popt[3])
        axs[i].plot(xx, func4(xx, *popt), '-',
                    label='H = %d, m=%.2e, c=%5.2f' % (h, popt[0], popt[1]))
    axs[i].legend()
    axs[i].set_ylabel('Log(Mass/Area)')
    axs[i].set_xlabel(r'Radius$^2$')
    axs[i].set_title(r'$\phi \in %s$' %
                     phi, fontsize=14)
    # axs[i].set_xlim(0, 1.7e10)
    # slopes.append(m)
    # intercepts.append(c)
plt.tight_layout()
plt.show()
# plt.savefig("./data/%s_trial_%d/line_fits_trial_%d.png" %
#             (param, trial, trial), dpi=200, format='png')
#
# for phi, slope in zip(phis, slopes):
#     plt.plot(vals, slope, label=phi, marker='o')
# plt.title(r'Plot of slope against column height', fontsize=14)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
# plt.xlabel("Column Height (km)")
# plt.ylabel("Slope")
# plt.tight_layout()
# plt.legend()
# # plt.show()
# plt.savefig("./data/%s_trial_%d/slope_trial_%d.png" %
#             (param, trial, trial), dpi=200, format='png')
#
# for phi, intercept in zip(phis, intercepts):
#     plt.plot(vals, intercept, label=r'\(\phi \in %s\)' %
#              phi, marker='o')
# plt.title(r'Plot of intercept against column height', fontsize=14)
# plt.xlabel("Column Height (km)")
# plt.ylabel("Intercept")
# plt.tight_layout()
# plt.legend()
# # plt.show()
# plt.savefig("./data/%s_trial_%d/intercept_trial_%d.png" %
#             (param, trial, trial), dpi=200, format='png')


fig, axs = plt.subplots(rows, 2, figsize=(
    10, 12), facecolor='w', edgecolor='k')
axs = axs.ravel()
slopes = []
intercepts = []
for i, phi in enumerate(phis):
    a = []
    b = []
    for h, erp in eruptions.items():
        print(str(h) + phi)
        new_df = erp.df
        ma = new_df['MassArea'].values
        wt = new_df[phi].values
        xx = new_df['radius'].values
        yy = ma * (wt / 100)

        pnts = list(set(zip(xx, yy)))

        pnts = sorted(pnts, key=lambda x: x[0])
        x = [p[0] for p in pnts if p[0] < 20000]
        y = [p[1] for p in pnts if p[0] < 20000]

        x = np.array(x)
        y = np.array(y)
        vals = np.polyfit(x, np.log(y), 1, w=np.log(y))
        vals
        yf = np.exp(vals[1]) * np.exp(x * vals[0])
        a.append(vals[0])
        b.append(vals[1])
        axs[i].plot(x, yf, '-',
                    label='H = %d, m=%.2e, c=%5.2f' % (h, vals[0], vals[1]))
    axs[i].legend()
    axs[i].set_ylabel('Log(Mass/Area)')
    axs[i].set_xlabel(r'Radius$^2$')
    axs[i].set_title(r'$\phi \in %s$' %
                     phi, fontsize=14)
    # axs[i].set_xlim(0, 1.7e10)
    # slopes.append(m)
    # intercepts.append(c)
plt.tight_layout()
plt.show()

new_df = eruptions[10000].df

ma = new_df['MassArea'].values
phi = new_df[phis[3]].values
xx = new_df['radius'].values
yy = ma * (phi / 100)

pnts = list(set(zip(xx, yy)))

pnts = sorted(pnts, key=lambda x: x[0])

ints = [np.trapz(yi, x=xi) for xi, yi in pnts]

x = np.array([p[0] for p in pnts])
y = np.array([p[1] for p in pnts])


vals = np.polyfit(x, np.log(y), 1, w=1 / x)
vals
yf = np.exp(vals[1]) * np.exp(x * vals[0])
yf
plt.plot(x, y, '.')
plt.plot(x, yf)
plt.savefig("polyfit_10000_y.png", dpi=200, format='png')

x = np.array(x)
y = np.array(y)
vals = np.polyfit(x, np.log(y), 1)
vals
yf = np.exp(vals[1]) * np.exp(x * vals[0])
yf
plt.plot(x, y, '.')
plt.plot(x, yf)
plt.show()
