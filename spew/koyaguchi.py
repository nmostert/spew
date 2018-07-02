from eruption import Eruption
from shapely.geometry import Point
import gridutils as grd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def radial_sample(eruption):
    # only works on regular grids for now

    ix_df = grd.get_index_grid(eruption.df)

    # Find grid index of pseudo origin.
    y_idx = ix_df.index.get_loc(eruption.grid_vent.y)
    x_idx = ix_df.columns.get_loc(eruption.grid_vent.x)

    # construct a radial line from the pseudo origin to the edge (East)
    line = eruption.df.iloc[ix_df.iloc[y_idx, x_idx:len(ix_df.columns)]]

    def dist(row):
        return eruption.vent_data.geometry.distance(row['geometry'])

    radius = line.apply(lambda row: dist(row), axis=1)
    line = line.assign(r=radius)
    return line


def plot_lnS_rsqr(samples, phi_labels):
    r_sqr = line['r'].values ** 2
    plt.xlabel("r^2")
    plt.ylabel("ln(S)")
    for phi in phi_labels:
        S = (samples['Mass/Area'].values * (samples[phi]))
        plt.plot(r_sqr, np.log(S), '.')
        plt.legend(phi_labels)


def plot_S_rsqr(samples, phi_labels):
    r_sqr = line['r'].values ** 2
    plt.xlabel("r^2")
    plt.ylabel("S")
    for phi in phi_labels:
        S = (line['Mass/Area'].values * (line[phi]))
        plt.plot(r_sqr, S, '.')
        plt.legend(phi_labels)


def coef_def(y, y_fit):
    # residual sum of squares
    ss_res = np.sum((y - y_fit) ** 2)

    # total sum of squares
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    # r-squared
    r2 = 1 - (ss_res / ss_tot)

    return r2


def minimise_coef_def(f, x, y, tol):
    fit = False
    r2 = 0

    while not fit:
        params = curve_fit(f, x, y, method='lm')
        yf = func(x, params[0][0], params[0][1])
        r2 = coef_def(y, yf)
        if r2 < tol and len(y) > 5:
            x = np.delete(x, 0)
            y = np.delete(y, 0)
            yf = np.delete(yf, 0)
        else:
            fit = True
    return x, y, yf, params


if __name__ == "__main__":
    filename = './data/cerroNegro_regGrid_noWind_SOURCE.txt'

    cn = Eruption(data=filename, vent=Point(532290, 1382690))

    # Find grid point closest to the vent (pseudo origin)

    line = radial_sample(cn)

    line['r_sqr'] = line.apply(lambda row: row.r ** 2, axis=1)

    def S_dict(row, phi_labels):
        S = dict([(phi, (row['Mass/Area'] * row[phi])) for phi in phi_labels])
        return S

    def logS_dict(row, phi_labels):
        logS = dict([(phi, np.log(row['Mass/Area'] * row[phi]))
                     for phi in phi_labels])
        return logS
    line['S'] = line.apply(lambda row: S_dict(row, cn.phi_labels), axis=1)

    line['logS'] = line.apply(
        lambda row: logS_dict(row, cn.phi_labels), axis=1)
    line.head()

    plot_lnS_rsqr(line, cn.phi_labels)
    plt.show()

    plot_S_rsqr(line, cn.phi_labels)
    plt.show()

    phi_df = pd.DataFrame(index=cn.phi_labels,
                          columns=['S', 'logS', 'clean_log_S', 'xx',
                                   'params', 'fit_vals'])

    phi_df.S = [line.apply(lambda row: row['S'][phi],
                           axis=1).values for phi in cn.phi_labels]

    def func(x, m, c):
        return m * x + c

    log_S_list = []
    clean_log_S_list = []
    xx_list = []
    params_list = []
    fit_vals_list = []
    r_sqr = line['r_sqr'].values
    for i, phi in enumerate(cn.phi_labels):
        ls = np.log(line['Mass/Area'].values * (line[phi]))
        ls = ls.replace([np.inf, -np.inf], np.nan)
        nan_idx = [not na for na in ls.isna()]
        cls = ls.dropna()
        xx = r_sqr[nan_idx]
        xx_list.append(xx)
        params = curve_fit(func, xx, cls, method='lm')
        fit_vals = func(xx, params[0][0], params[0][1])

        lin_x, lin_cls, lin_fit_vals, lin_params = minimise_coef_def(
            func, xx, cls.values, 0.9999)
        residuals = cls - fit_vals
        fit_res = lin_cls - lin_fit_vals
        # plt.xlabel("r^2")
        # plt.ylabel("ln(S)")
        f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        ax1.scatter(xx, cls, marker='.')
        ax1.plot(xx, fit_vals, linestyle='-', color='orange')
        ax1.plot(lin_x, lin_fit_vals, linestyle=':', color='green')
        lab = 'y = ' + '{:0.3e}'.format(params[0][0]) + \
            'x + ' + str(round(params[0][1], 4))
        lin_lab = 'y = ' + '{:0.3e}'.format(lin_params[0][0]) + \
            'x + ' + str(round(lin_params[0][1], 4))
        ax1.legend([lab, lin_lab, phi])
        ax1.set_title('Linear fit')
        ax1.set_xlabel('r^2')
        ax1.set_ylabel('ln(S)')

        print(coef_def(lin_cls, lin_fit_vals))

        ax2.vlines(xx, [0], residuals,
                   color='orange', linestyle='-')
        ax2.vlines(lin_x, [0], fit_res, color='g', linestyle=':')
        ax2.set_xlabel('r^2')
        ax2.set_ylabel('Residual')
        ax2.set_title('Residuals')
        red_lab = "R^2 = " + str(round(coef_def(cls.values, fit_vals), 6))
        lin_red_lab = "R^2 = " + str(round(coef_def(lin_cls, lin_fit_vals), 6))
        ax2.legend([red_lab, lin_red_lab])
        plt.savefig('res_' + str(i) + '.eps', dpi=200, format='eps')

        log_S_list.append(ls.values)
        clean_log_S_list.append(cls.values)
        params_list.append(params)
        fit_vals_list.append(fit_vals)

    phi_df.logS = log_S_list
    phi_df.clean_log_S = clean_log_S_list
    phi_df.xx = xx_list
    phi_df.params = params_list
    phi_df.fit_vals = fit_vals_list

    phi_df
