"""
log and exp used for evaluating formulas
errorAnalysis for error analysis
arange and meshgrid for 3d plot coordinates
matplotlib.pyplot for plotting
Axes3D for 3d axes
cm for colormaps
"""
from math import log, exp
import errorAnalysis as err
from numpy import arange, meshgrid
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

"""This python modual containts functions for the calculation of thermophysical properties of moist
air in the range between 0 and 100 degC. For details regarding these functions see P.T. Tsilingiris
(2008) "Thermophsical and transport properties of humid air at temperature range between 0 and 100
degC," Energy Conv. Manag. (49) pp. 1098-1110.

Author: Rich H. Inman
email: rinman@ucsd.edu
Date Created: 21 Nov. 2016
"""

def build_xy_grid(x_det, y_det):
    """Use x and y grid details to return meshgrid

    This function builds a meshgrid to evaluate and plot functions of two variables.

    Positional argument(s)
    x_det -- [x_0, x_f, delta_x]
    y_det -- [y_0, y_f, delta_y]

    Output(s):
    xgrid, ygrid -- lists of lists for (x,y) coordinates
    """
    xgrid = arange(x_det[0], x_det[1]+x_det[2], x_det[2])
    ygrid = arange(y_det[0], y_det[1]+y_det[2], y_det[2])
    return meshgrid(xgrid, ygrid)

def plot_property(x_grid, y_grid, r_hum=0.5, ver=0):
    """Use x_grid and y_grid to plot a property(T, P, RH).
    
    Preliminary function for plotting thermophysical properties of moist air. Much of this code will
    end up being recycled and wrapped up into neater functions.

    Positional arguments(s)
    x_grid -- x-numpy.meshgrid
    y_grid -- y-numpy.meshgrid

    Output(s):
    None
    """
    z_grid = x_grid*float('nan')
    fig = plt.figure()
    axis = fig.gca(projection='3d')
    verts = []
    for i, _ in enumerate(x_grid):
        for j, _ in enumerate(x_grid[0]):
            if ver == 0: # p_sat
                z_grid[i][j] = eff_p_sat_liq(x_grid[i][j], y_grid[i][j])
            elif ver == 1: # f_en
                z_grid[i][j] = enh_fact_liq(x_grid[i][j], y_grid[i][j], p_sat_liq(x_grid[i][j]))
            elif ver == 2: # x_v
                z_grid[i][j] = r_hum*eff_p_sat_liq(x_grid[i][j], y_grid[i][j])/y_grid[i][j]
                if z_grid[i][j] > 1:
                    z_grid[i][j] = 1
    axis.plot_surface(x_grid, y_grid, z_grid, rstride=4, cstride=4, alpha=0.3)
    axis.contour(x_grid, y_grid, z_grid, zdir='x', offset=366, cmap=cm.coolwarm)
    axis.contour(x_grid, y_grid, z_grid, zdir='y', offset=112500, cmap=cm.coolwarm)
    axis.set_xlabel('\n'+r'$T$, [K]', linespacing=1.8)
    axis.set_xticks([260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360])
    axis.set_xticklabels(['', 270, '', 290, '', 310, '', 330, '', 350, ''])
    axis.set_xlim(366,260)
    axis.set_ylabel('\n'+r'$P$, [hPa]', linespacing=2)
    axis.set_yticks([25000, 37500, 50000, 62500, 75000, 87500, 100000])
    axis.set_yticklabels([250, '', 500, '', 750, '', 1000])
    axis.set_ylim(12500, 112500)
    if ver == 0: # p_sat
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=0, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'$P_{\rm{sat}}$, [hPa]', linespacing=1.5)
        axis.set_zticks([0, 10000, 20000, 30000, 40000, 50000])
        axis.set_zticklabels([0, 100, 200, 300, 400, 500])
        axis.set_zlim(0, 50000)
        plt.savefig('000_P_sat(T,P).pdf')
    elif ver == 1: # f_en
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=0.992, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'f$_{\rm{en}}$', linespacing=2.5, style='italic')
        axis.tick_params(axis='z', which='major', pad=8)
        axis.set_zlim(0.992, 1.006)
        plt.savefig('001_f_en(T,P).pdf')
    elif ver == 2: # x_v
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=0, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'x$_v$', linespacing=1.5, style='italic')
        axis.set_zticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        axis.set_zlim(0, 1)
        axis.text(250, 0, 1.25, 'RH = '+str(int(r_hum*100))+'%')
        plt.savefig('002_x_v(T,P)'+str(int(r_hum*100))+'.pdf')

    plt.show()

def get_temp_k_obs(temp_k, delta_t=0.2):
    """Use temp k to return list of tuples containing uncertainties.

    This function accepts either a single observation or an iterable of observations and appends the
    uncertainty to each of the observations. The default uncertainty is 0.2 K.

    Positional argument(s):
    temp_k  -- temperature(s) in Kelvin

    Keyword argument(s):
    delta_t -- uncertainty of all temperature observations (default 0.2)

    Output(s):
    temp_k_obs -- list of tuples [(obs_1, delta_1), ..., (obs_n, delta_n)]
    """
    try:
        temp_k_obs = [(float(temp), delta_t) for temp in temp_k]
    except TypeError:
        temp_k_obs = [(float(temp_k), delta_t)]
    return temp_k_obs

def get_p_pa_obs(p_pa, frac_p=0.0025):
    """Use p_pa to return list of tuples containing uncertainties.

    This function accepts either a single observation or an iterable of observations and appends the
    uncertainty to each of the observations. The default fractional uncertainty is 0.25% or 0.0025.

    Positional argument(s):
    p_pa  -- pressure(s) in Pascal

    Keyword argument(s):
    frac_p -- fractional uncertainty of all pressure observations (default 0.25% or 0.0025)

    Output(s)
    p_pa_obs -- list of tuples [(obs_1, delta_1), ..., (obs_n, delta_n)]
    """
    try:
        p_pa_obs = [(float(p), frac_p*p) for p in p_pa]
    except TypeError:
        p_pa_obs = [(float(p_pa), frac_p*p_pa)]
    return p_pa_obs

def check_temp_liq(temp_k):
    """Check the temperature in Kelvin to return boolean of liquid water

    This funciton checkes the temperature in Kelvin to see if the temperature is between 273.15 and
    373.15 K. Returns True if the temperature is within said range and False otherwise. This
    function assumes that the standard form of an observation has been passed to it. That is, a
    list of tuples containing the best guess at the observation in the [0] index and the uncertainty
    of the observation in the [1] index.

    Positional argument(s):
    temp_k -- list of tuples [(temp_k, delta_t)] in Kelvin
    
    Output(s):
    if temperature is between 273.15 and 373.15 K:
        True
    else:
        False
    """
    if temp_k[0] < 273.15:
        print '\tError in p_sat_l: Temperature below 273.15 K.'
        return False
    elif temp_k[0] > 373.15 + temp_k[1]:
        print '\tError in p_sat_l: Temperature above 373.15 K.'
        return False
    else:
        return True

def p_sat_liq(temp_k):
    """Use Hardy (1998) to return P_sat(T_K) in Pa.

    This function uses Hardy (1998), "ITS-90 Formulations for Vapor Pressure, Frostpoint
    Temperature, Dewpoint Temperature, and Enhancement Factors in the Range -100 to +100 C," to
    calculate the saturation vapor pressure over liquid water in the range 0 to 100 C.

    Positional argument(s):
    temp_k -- temperature in Kelvin
    
    Output(s):
    p_sat  -- saturation vapor pressure over liquid water in Pa
    """
    g_coef = (-2.8365744e3, -6.028076559e3, 1.954263612e1, -2.737830188e-2, 1.6261698e-5,
              7.0229056e-10, -1.8680009e-13, 2.7150305)
    p_sat = 0
    for i in range(7):
        p_sat += g_coef[i] * pow(temp_k, i-2)
    p_sat += g_coef[7] * log(temp_k)
    p_sat = exp(p_sat)
    return p_sat

def enh_fact_liq(temp_k, p_pa, p_sat):
    """Use Hardy (1998) to return f_enh.

    This function uses Hardy (1998), "ITS-90 Formulations for Vapor Pressure, Frostpoint
    Temperature, Dewpoint Temperature, and Enhancement Factors in the Range -100 to +100 C," to
    calculate the enhancement factor over liquid water in the range 0 to 100 C.

    Positional argument(s):
    temp_k -- temperature in Kelvin
    p_pa   -- total pressure in Pa
    p_sat  -- saturation vapor pressure of water over liquid water in Pa
    
    Output(s):
    f_enh -- enhancement factor over liquid water
    """
    a_coef = (-1.6302041e-1, 1.8071570e-3, -6.7703064e-6, 8.5813609e-9)
    alpha = 0
    b_coef = (-5.9890467e1, 3.4378043e-1, -7.7326396e-4, 6.3405286e-7)
    beta = 0
    for i in range(4):
        alpha += a_coef[i] * pow(temp_k, i)
        beta += b_coef[i] * pow(temp_k, i)
    beta = exp(beta)
    f_enh = exp(alpha*(1 - p_sat/p_pa) + beta*(p_pa/p_sat - 1))
    return f_enh

def eff_p_sat_liq(temp_k, p_pa):
    """Use p_sat_liq and enh_fact_liq to return effective saturation pressure.

    This function uses Hardy (1998), "ITS-90 Formulations for Vapor Pressure, Frostpoint
    Temperature, Dewpoint Temperature, and Enhancement Factors in the Range -100 to +100 C," to
    calculate the effective saturation vapor pressure over liquid water in the range 0 to 100 C.

    Positional argument(s):
    temp_k -- temperature in Kelvin
    p_pa   -- total pressure in Pa
    
    Output(s):
    eff_p  -- effective saturation vapor pressure over liquid water
    """
    p_sat = p_sat_liq(temp_k)
    enh_fact = enh_fact_liq(temp_k, p_pa, p_sat)
    return enh_fact*p_sat

def main():
    """Main"""