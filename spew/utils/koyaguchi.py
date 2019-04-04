from ..core.eruption import Eruption, get_index_grid
from shapely.geometry import Point
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import comb


def radial_sample(eruption):
    """

    Parameters
    ----------
    eruption: Spew Eruption
        Eruption class instance

    Returns
    -------
    line: Pandas Dataframe
        Dataframe containing points along a line leading directly East from vent.

    """
    ix_df = get_index_grid(eruption.df)

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
    plt.xlabel(r"$r^2$")
    plt.ylabel(r"$\ln(S)$")
    for phi in phi_labels:
        S = (samples['MassArea'].values * (samples[phi]))
        plt.plot(r_sqr, np.log(S), '.')
        plt.legend(phi_labels)


def plot_S_rsqr(samples, phi_labels):
    r_sqr = line['r'].values ** 2
    plt.xlabel(r"$r^2$")
    plt.ylabel(r"$S$")
    for phi in phi_labels:
        S = (line['MassArea'].values * (line[phi]))
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


def proximal_deletion(f, x, y, tol):
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


def phi_to_d(phi):
    return 0.001 * (2 ** -phi)


def d_to_phi(d):
    return -np.log2(d / 0.001)


def high_re_velocity(phi, clast_density, air_density, drag_coeff=1):
    g = 9.80655
    v_h = drag_coeff * \
        np.sqrt((clast_density * g * phi_to_d(phi)) / air_density)
    return v_h


def low_re_velocity(phi, clast_density, air_density, air_viscosity):
    g = 9.80655
    v_l = ((clast_density - air_density) * g *
           (phi_to_d(phi)**2)) / (18 * air_viscosity)
    return v_l


def re_intercept_grainsize(clast_density, air_density, air_viscosity, drag_coeff=1):
    """Phi intercept between low and high Re settling velocities.

    Grain-size at which the low Re settling velocity formula intersects the
    high Re settling velocity formula.
    """
    g = 9.80655
    d = ((18 * air_viscosity * drag_coeff * np.sqrt((clast_density * g) /
                                                    air_density)) / ((clast_density - air_density) * g)) ** (2 / 3)
    return d_to_phi(d)


def tanh_smoother(phi, A, B):
    s = 0.5 + 0.5 * np.tanh((phi - A) / B)
    return s


def smoothed_function(f, g, s):
    h = s * g + (1 - s) * f
    return h


def bernstein_poly(i, n, t):
    """
    The Bernstein polynomial of n, i as a function of t
    """

    return comb(n, i) * (t**(n - i)) * (1 - t)**i


def bezier_curve(points, nTimes=1000):
    """
       Given a set of control points, return the
       bezier curve defined by the control points.

       points should be a list of lists, or list of tuples
       such as [ [1,1],
                 [2,3],
                 [4,5], ..[Xn, Yn] ]
        nTimes is the number of time steps, defaults to 1000

        See http://processingjs.nihongoresources.com/bezierinfo/
    """

    nPoints = len(points)
    xPoints = np.array([p[0] for p in points])
    yPoints = np.array([p[1] for p in points])

    t = np.linspace(0.0, 1.0, nTimes)

    polynomial_array = np.array(
        [bernstein_poly(i, nPoints - 1, t) for i in range(0, nPoints)])

    xvals = np.dot(xPoints, polynomial_array)
    yvals = np.dot(yPoints, polynomial_array)

    return xvals, yvals


if __name__ == "__main__":
    filename1 = './data/ceroNegro_GITHUB.txt'
    cn = Eruption(data=filename1, vent=Point(532290, 1382690), test=False)

    # Find grid point closest to the vent (pseudo origin)

    line = radial_sample(cn)

    line['r_sqr'] = line.apply(lambda row: row.r ** 2, axis=1)

    def S_dict(row, phi_labels):
        S = dict([(phi, (row['MassArea'] * row[phi])) for phi in phi_labels])
        return S

    def logS_dict(row, phi_labels):
        logS = dict([(phi, np.log(row['MassArea'] * row[phi]))
                     for phi in phi_labels])
        return logS
    line['S'] = line.apply(lambda row: S_dict(row, cn.phi_labels), axis=1)

    line['logS'] = line.apply(
        lambda row: logS_dict(row, cn.phi_labels), axis=1)
    line.head()

    plot_lnS_rsqr(line, cn.phi_labels)
    # plt.savefig('lnS.eps', dpi=200, format='eps')
    # plt.show()

    plot_S_rsqr(line, cn.phi_labels)
    # plt.savefig('res_' + str(i) + '.eps', dpi=200, format='eps')
    # plt.show()

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
        ls = np.log(line['MassArea'].values * (line[phi]))
        ls = ls.replace([np.inf, -np.inf], np.nan)
        nan_idx = [not na for na in ls.isna()]
        cls = ls.dropna()
        xx = r_sqr[nan_idx]
        xx_list.append(xx)
        params = curve_fit(func, xx, cls, method='lm')
        fit_vals = func(xx, params[0][0], params[0][1])

        lin_x, lin_cls, lin_fit_vals, lin_params = proximal_deletion(
            func, xx, cls.values, 0.9999)
        residuals = cls - fit_vals
        fit_res = lin_cls - lin_fit_vals
        # plt.xlabel("r^2")
        # plt.ylabel("ln(S)")
        f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        ax1.scatter(xx, cls, marker='.')
        ax1.plot(xx, fit_vals, linestyle='-', color='orange')
        ax1.plot(lin_x, lin_fit_vals, linestyle=':', color='green')
        lab = r'$y = $' + '{:0.3e}'.format(params[0][0]) + \
            r'$x + $' + str(round(params[0][1], 4))
        lin_lab = r'$y = $' + '{:0.3e}'.format(lin_params[0][0]) + \
            r'$x + $' + str(round(lin_params[0][1], 4))
        ax1.legend([lab, lin_lab, phi])
        ax1.set_title('Linear fit')
        ax1.set_xlabel(r'$r^2$')
        ax1.set_ylabel(r'$\ln(S)$')

        print(len(xx), len(lin_x))

        ax2.vlines(xx, [0], residuals,
                   color='orange', linestyle='-')
        ax2.vlines(lin_x, [0], fit_res, color='g', linestyle=':')
        ax2.set_xlabel(r'$r^2$')
        ax2.set_ylabel('Residual')
        ax2.set_title('Residuals')
        red_lab = r"$R^2 = $" + str(round(coef_def(cls.values, fit_vals), 6))
        lin_red_lab = r"$R^2 = $" + \
            str(round(coef_def(lin_cls, lin_fit_vals), 6))
        ax2.legend([red_lab, lin_red_lab])
        # plt.savefig('res_' + str(i) + '.eps', dpi=200, format='eps')

        log_S_list.append(ls.values)
        clean_log_S_list.append(cls.values)
        params_list.append(params)
        fit_vals_list.append(fit_vals)

    phi_df.logS = log_S_list
    phi_df.clean_log_S = clean_log_S_list
    phi_df.xx = xx_list
    phi_df.params = params_list
    phi_df.fit_vals = fit_vals_list

    phi_centers = np.array([-4, -3.5, -3, -2.5, -2, -1.5, -
                            1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4])
    clast_density = np.array([1000, 1000, 1000, 1200, 1100, 1100, 1200, 1500,
                              1600, 1800, 2000, 2200, 2300, 2600, 2900, 2400,
                              2400])

    # clast_density = np.array([2700] * 17)

    # things done at 25 km altitude
    air_density = 3.69 * (10**-2)
    air_viscosity = 1.40 * (10**-5)
    args_high = [air_density, 1]
    args_low = [air_density, air_viscosity]

    f = low_re_velocity(phi_centers, clast_density, *args_low)
    g = high_re_velocity(phi_centers, clast_density, *args_high)

    A = re_intercept_grainsize(clast_density, air_density, air_viscosity)

    s = tanh_smoother(phi_centers, A, 1)

    h = smoothed_function(f, g, s)

    low_ints = [low_re_velocity(ints, dens, *args_low)
                for ints, dens in zip(A, clast_density)]
    high_ints = [high_re_velocity(ints, dens, *args_high)
                 for ints, dens in zip(A, clast_density)]

    low_ints
    high_ints
    plt.figure()
    plt.semilogy(phi_centers, f, 'k--')
    plt.plot(phi_centers, g, 'g--')
    ax = plt.gca()
    ax.invert_xaxis()
    plt.xlabel(r"Grain size ($\Phi$-scale)")
    plt.ylabel("Terminal velocity (m/s)")
    # plt.xlim([0,3])
    # plt.ylim([10,100])
    plt.title("Terminal velocity in Low and High Reynolds domains")
    plt.legend(["Low Re", "High Re", "Smoothed function"])
    plt.savefig('term_vel.eps', dpi=200, format='eps')
    plt.show()
    A
    print(low_ints)

    phi_centers

    lf = np.log(f)
    lg = np.log(g)

    points = [(phi_centers[0], lg[0]), (phi_centers[10], lf[10]), (phi_centers[11], lg[11]),
              (phi_centers[16], lf[16])]

    points
    xpoints = [p[0] for p in points]
    ypoints = [p[1] for p in points]
    xvals, yvals = bezier_curve(points, nTimes=1000)

    plt.figure()
    plt.semilogy(phi_centers, f, 'k--')
    plt.plot(phi_centers, g, 'g--')
    plt.plot(xvals, np.exp(yvals))
    plt.plot(xpoints, np.exp(ypoints), "ro")
    plt.xlabel(r"Grain size ($\Phi$-scale)")
    plt.ylabel("Terminal velocity (m/s)")
    plt.legend(["Low Re", "High Re", r"B\'ezier curve", "Weights"])
    plt.title(r"B\'ezier curve transition between Low and High Reynolds domains")
    for nr in range(len(points)):
        plt.text(points[nr][0], np.exp(points[nr][1]), nr)
    ax = plt.gca()
    ax.invert_xaxis()
    plt.savefig('bez_vel.eps', dpi=200, format='eps')
    plt.show()
