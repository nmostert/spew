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

###################OPTION ONE: WEIGHTED POLYFIT##########################

for k, phi in enumerate(phis):
    fig, axs = plt.subplots(rows, 2, figsize=(
        10, 15), facecolor='w', edgecolor='k')
    axs = axs.ravel()
    for i, (h, erp) in enumerate(eruptions.items()):
        new_df = erp.df
        ma = new_df['MassArea'].values
        wt = new_df[phi].values
        xx = new_df['radius'].values
        yy = ma * (wt / 100)

        pnts = list(set(zip(xx, yy)))
        pnts = sorted(pnts, key=lambda x: x[0])

        x = np.array([p[0] for p in pnts])
        y = np.array([p[1] for p in pnts])

        tot_mass = np.trapz(y, x=x)
        idx = 0
        perc = 0.0
        while perc < 99.99999999:
            mass = np.trapz(y[:idx], x=x[:idx])
            idx += 1
            perc = (mass / tot_mass) * 100

        x = x[:idx]
        y = y[:idx]

        vals_1 = np.polyfit(x, np.log(y), 1, w=y)
        vals_2 = np.polyfit(x, np.log(y), 1, w=1 / x)
        vals_3 = np.polyfit(x, np.log(y), 1, w=1 / (x**2))
        xf = np.linspace(x.min(), x.max(), 100)
        yf_1 = np.exp(vals_1[1]) * np.exp(xf * vals_1[0])
        yf_2 = np.exp(vals_2[1]) * np.exp(xf * vals_2[0])
        yf_3 = np.exp(vals_3[1]) * np.exp(xf * vals_3[0])
        axs[i].scatter(x, y, s=30, facecolors='none', edgecolors='r')
        axs[i].plot(xf, yf_1, ':', color='k', lw=1.2, label=r'$w=y$, $a$=%.2e, $b$=%.2e' % (
            np.exp(vals_1[1]), vals_1[0]))
        axs[i].plot(xf, yf_2, '--', color='k', lw=1.2, label=r'$w=1/x$, $a$=%.2e, $b$=%.2e' %
                    (np.exp(vals_2[1]), vals_2[0]))
        axs[i].plot(xf, yf_3, '-.', color='k', lw=1.2, label=r'$w=1/x^2$, $a$=%.2e, $b$=%.2e' %
                    (np.exp(vals_3[1]), vals_3[0]))
        axs[i].legend()
        axs[i].set_ylabel('Mass/Area')
        axs[i].set_xlabel(r'Radius')
        axs[i].set_title(r'$\phi \in %s$, $H = %d$km' %
                         (phi, disp_func(h)), fontsize=14)
    plt.tight_layout()
    plt.savefig("./data/%s_trial_%d/polyfit_weights_%d.pdf" %
                (param, trial, k), format='pdf')


###################SCIPY.MINIMIZE##########################


def gamma_kernel(x, a, b, c):
    return c * (x**a) * e**(b * x)


def sse(k, x_data, y_data):
    return sum((y_data - gamma_kernel(x_data, *k)) ** 2)


fit_df = pd.DataFrame(columns=['H', 'Phi', 'a', 'b', 'c'])

for k, phi in enumerate(phis):
    fig, axs = plt.subplots(rows, 2, figsize=(
        10, 15), facecolor='w', edgecolor='k')
    axs = axs.ravel()
    for i, (h, erp) in enumerate(eruptions.items()):
        new_df = erp.df
        ma = new_df['MassArea'].values
        wt = new_df[phi].values
        xx = new_df['radius'].values
        yy = ma * (wt / 100)

        pnts = list(set(zip(xx, yy)))
        pnts = sorted(pnts, key=lambda x: x[0])

        x = np.array([p[0] for p in pnts])
        y = np.array([p[1] for p in pnts])

        tot_mass = np.trapz(y, x=x)
        idx = 0
        perc = 0.0
        while perc < 99.99999999:
            mass = np.trapz(y[:idx], x=x[:idx])
            idx += 1
            perc = (mass / tot_mass) * 100

        x_data = x[:idx]
        y_data = y[:idx]

        k0 = (0, 0, 0)

        def fun(k): return gamma_kernel_sse(k, x_data, y_data)

        vals = minimize(fun, k0)

        xf = np.linspace(x_data.min(), x_data.max(), 100)
        yf = gamma_kernel(xf, *vals.x)

        axs[i].scatter(x_data, y_data, s=30, facecolors='none', edgecolors='r')
        axs[i].plot(xf, yf, '--', color='k', lw=1.2, label=r'\noindent $a$=%.2e,\\ $b$=%.2e,\\ $c$=%.2e' % (
            vals.x[0], vals.x[1], vals.x[2]))
        axs[i].legend()
        axs[i].set_ylabel('Mass/Area')
        axs[i].set_xlabel(r'Radius')
        axs[i].set_title(r'$\phi \in %s$, $H = %d$km' %
                         (phi, disp_func(h)), fontsize=14)
        fit_df = fit_df.append(
            {'H': h, 'Phi': phi, 'a': vals.x[0], 'b': vals.x[1], 'c': vals.x[2]}, ignore_index=True)
    plt.tight_layout()
    plt.show()
    # plt.savefig("./data/%s_trial_%d/minimize_%d.pdf" %
    #             (param, trial, k), format='pdf')

fit_df = fit_df[fit_df['H'] != 5000]

fig, axs = plt.subplots(7, 3, figsize=(
    10, 15), facecolor='w', edgecolor='k', sharex=True)
axs = axs.ravel()
axi = 0
H = list(eruptions.keys())
for h in H[1:]:
    group = fit_df[fit_df['H'] == h].groupby("Phi").mean().reindex(phis)

    for p, col, format in zip(['a', 'b', 'c'], ['steelblue', 'indianred', 'dimgrey'], ['%.2f', '%3f', '%.4f']):
        group[p].plot(kind="line", color=col, ax=axs[axi], legend=False)
        axs[axi].get_yaxis().get_major_formatter().set_powerlimits((-3, 4))
        axs[axi].set_xlabel(r'Grainsize ($\phi$)')
        if p == 'c':
            axs[axi].set_ylabel(r'Height = %d km' % disp_func(h))
            axs[axi].get_yaxis().set_label_position("right")
        if axi < 3:
            axs[axi].set_title(p)
        if axi < 20:
            axi += 1
plt.tight_layout()
plt.show()
# plt.savefig("./data/%s_trial_%d/h_by_phi_abc.pdf" %
#             (param, trial), format='pdf')

fig, axs = plt.subplots(8, 3, figsize=(
    10, 15), facecolor='w', edgecolor='k', sharex=True)
axs = axs.ravel()
len(axs)
axi = 0
for phi in phis:
    group = fit_df[fit_df['Phi'] == phi].groupby("H").mean()
    group.index = [10, 15, 20, 25, 30, 35, 40]
    for p, col in zip(['a', 'b', 'c'], ['steelblue', 'indianred', 'dimgrey']):
        group[p].plot(kind="bar", color=col, ax=axs[axi], legend=False)
        axs[axi].get_yaxis().get_major_formatter().set_powerlimits((-2, 3))
        axs[axi].set_xlabel(r'Column Height (km)')
        if p == 'c':
            axs[axi].set_ylabel(r'$\phi \in %s$' % phi)
            axs[axi].get_yaxis().set_label_position("right")
        if axi < 3:
            axs[axi].set_title(p)
        if axi < 23:
            axi += 1
plt.tight_layout()
plt.show()
# plt.savefig("./data/%s_trial_%d/phi_by_h_abc.pdf" %
#             (param, trial), format='pdf')

fit_df

a_piv = fit_df.pivot(index='Phi', columns='H', values='a')
b_piv = fit_df.pivot(index='Phi', columns='H', values='b')
c_piv = fit_df.pivot(index='Phi', columns='H', values='c')
a_piv = a_piv.reindex(phis)
b_piv = b_piv.reindex(phis)
c_piv = c_piv.reindex(phis)

xlabels = [123, 5, 10, 15, 20, 25, 30, 35, 40]
ylabels = ['213'] + phis


fig, axs = plt.subplots(3, 1, figsize=(
    10, 15), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, piv, p in zip([0, 1, 2], [a_piv, b_piv, c_piv], ['a', 'b', 'c']):
    c = axs[i].imshow(scipy.ndimage.zoom(piv.values, 3),
                      cmap="viridis", interpolation="spline16")
    cont = axs[i].contour(scipy.ndimage.zoom(piv.values, 3), colors='w')
    axs[i].set_xticklabels(xlabels)
    axs[i].set_yticklabels(ylabels)
    fig.colorbar(c, ax=axs[i])
    axs[i].grid(False)
    axs[i].set_title(p)
    axs[i].clabel(cont, inline=1, fontsize=10)
    axs[i].set_xlabel("Column height (km)")
    axs[i].set_ylabel("Phi class")
plt.savefig("./data/%s_trial_%d/param_maps_cont.pdf" %
            (param, trial), format='pdf')

piv

fig, axs = plt.subplots(3, 1, figsize=(
    10, 15), facecolor='w', edgecolor='k', subplot_kw=dict(projection='3d'))
axs = axs.ravel()
for i, piv, p in zip([0, 1, 2], [a_piv, b_piv, c_piv], ['a', 'b', 'c']):
    c = axs[i].plot_surface(range(8), piv.columns.astype(
        'float'), piv.values.astype('float'), cmap="hot")
    axs[i].set_xticklabels(xlabels)
    axs[i].set_yticklabels(ylabels)
    fig.colorbar(c, ax=axs[i])
    axs[i].grid(False)
    axs[i].set_title(p)
    axs[i].clabel(cont, inline=1, fontsize=10)
    axs[i].set_xlabel("Column height (km)")
    axs[i].set_ylabel("Phi class")
plt.show()
# plt.savefig("./data/%s_trial_%d/param_maps_cont.pdf" %
#             (param, trial), format='pdf')
