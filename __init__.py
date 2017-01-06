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
from matplotlib import cm, colors

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

def shave_xy_grid(x_grid, y_grid):
    """Use x_grid and y_grid to shave meshgrid

    This function sets coordinate paris in x_grid and y_grid to NaN if they satisfy the condition:
    P < 101325*pow(T/313.15, 10).

    Positional argument(s)
    x_grid -- meshgrid from build_xy_grid function
    y_grid -- meshgrid from build_xy_grid function

    Output(s):
    xgrid, ygrid -- lists of lists for (x,y) coordinates
    """
    for i, _ in enumerate(x_grid):
        for j, _ in enumerate(x_grid[0]):
            if y_grid[i][j] < 101325.*pow(x_grid[i][j]/313.15, 10):
                x_grid[i][j] = float('nan')
                y_grid[i][j] = float('nan')
    return x_grid, y_grid

def plot_property(x_grid, y_grid, r_hum=0.5, ver=0):
    """Use x_grid and y_grid to plot a property(T, P, RH).

    Preliminary function for plotting thermophysical properties of moist air. Much of this code will
    end up being recycled and wrapped up into neater functions.

    Positional arguments(s)
    x_grid -- x-numpy.meshgrid
    y_grid -- y-numpy.meshgrid

    Keyword arguments(s):
    r_hum  -- relative humidity as a fraction
    ver    -- version of the plot to produce
              {0: p_sat, 1: f_en, 2: x_v, 3: comp_fact, 4: rho_mix}

    Output(s):
    None
    """

    # Normalize colorbar
    if ver == 0: # p_sat
        lev = arange(0, 8100, 81)
    elif ver == 1: # f_en
        lev = arange(1.001, 1.00501, 5e-5)
    elif ver == 2: #x_v
        lev = arange(0, 0.081, 8.1e-4)
    elif ver == 3: #comp_fact
        lev = arange(0.997, 1, 3e-5)
    elif ver == 4: # rho_m
        lev = arange(0.3, 1.3, 0.001)
    elif ver == 5: # m_v
        lev = arange(0, 0.051, 5e-4)
    elif ver == 6: # rho_v
        lev = arange(0, 0.051, 5e-4)
    norml = colors.BoundaryNorm(lev, 256)

    # Initalize z-grid
    z_grid = x_grid*float('nan')

    # Create 3d figure
    fig = plt.figure(figsize=[9, 6])
    axis = fig.gca(projection='3d')

    # Build plot
    for i, _ in enumerate(x_grid):
        for j, _ in enumerate(x_grid[0]):
            if ver == 0: # p_sat
                z_grid[i][j] = eff_p_sat_liq(x_grid[i][j], y_grid[i][j])
            elif ver == 1: # f_en
                z_grid[i][j] = enh_fact_liq(x_grid[i][j], y_grid[i][j], p_sat_liq(x_grid[i][j]))
            elif ver == 2: # x_v
                z_grid[i][j] = x_v(x_grid[i][j], y_grid[i][j], r_hum)
            elif ver == 3: # comp_fact
                z_grid[i][j] = comp_fact(x_grid[i][j], y_grid[i][j])
            elif ver == 4: # rho_m
                z_grid[i][j] = rho_mix(x_grid[i][j], y_grid[i][j], r_hum)
            elif ver == 5: # m_v
                z_grid[i][j] = m_v(x_grid[i][j], y_grid[i][j], r_hum)
            elif ver == 6: # rho_v
                z_grid[i][j] = rho_v(x_grid[i][j], y_grid[i][j], r_hum)

    # Plot
    surf = axis.plot_surface(x_grid, y_grid, z_grid, rstride=2, cstride=2, cmap=cm.coolwarm,
                             linewidth=0, norm=norml, rasterized=True)
    axis.contour(x_grid, y_grid, z_grid, zdir='x', offset=320, cmap=cm.coolwarm)
    axis.contour(x_grid, y_grid, z_grid, zdir='y', offset=114025, cmap=cm.coolwarm)

    # x-axis
    axis.set_xlabel('\n'+r'$T$, [K]', linespacing=1.8)
    axis.set_xticks([270, 280, 290, 300, 310])
    axis.set_xlim(320, 266)

    # y-axis
    axis.set_ylabel('\n'+r'$P$, [hPa]', linespacing=2)
    axis.set_yticks([20e3, 40e3, 60e3, 80e3, 100e3])
    axis.set_yticklabels([200, 400, 600, 800, '1,000'])
    axis.set_ylim(13200, 114025)

    # z-axis
    if ver == 0: # p_sat
        z_ticks = [0, 2000, 4000, 6000, 8000]
        fig.colorbar(surf, shrink=0.5, aspect=10, pad=0.07, format='%.0f', ticks=z_ticks)
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=-523, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'$P_{\rm{sat}}$, [hPa]', linespacing=1.5)
        axis.set_zticks(z_ticks)
        axis.set_zticklabels([0, 20, 40, 60, 80])
        axis.set_zlim(-523, 8555)
        plt.savefig('000_P_sat(T,P).pdf')
    elif ver == 1: # f_en
        z_ticks = [1.001, 1.002, 1.003, 1.004, 1.005]
        fig.colorbar(surf, shrink=0.5, aspect=10, pad=0.07, format='%.3f', ticks=z_ticks)
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=1.0006, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'f$_{\rm{en}}$', linespacing=3.2, style='italic')
        axis.tick_params(axis='z', which='major', pad=8)
        axis.set_zticks(z_ticks)
        axis.set_zticklabels(z_ticks)
        axis.set_zlim(1.0006, 1.0053)
        plt.savefig('001_f_en(T,P).pdf')
    elif ver == 2: # x_v
        z_ticks = [0, 0.02, 0.04, 0.06, 0.08]
        fig.colorbar(surf, shrink=0.5, aspect=10, pad=0.07, format='%.2f', ticks=z_ticks)
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=-0.01, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'x$_v$', style='italic')
        axis.set_zticks(z_ticks)
        axis.set_zlim(-0.01, 0.085)
        axis.text(285, 100000, 0.06, 'RH = '+str(int(r_hum*100))+'%')
        plt.savefig('002_x_v(T,P)'+str(int(r_hum*100))+'.pdf')
    elif ver == 3: # comp_fact
        z_ticks = [0.997, 0.9975, 0.998, 0.9985, 0.999, 0.9995, 1]
        fig.colorbar(surf, shrink=0.5, aspect=10, pad=0.07, format='%.3f', ticks=z_ticks)
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=0.997, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'Z$_c$', style='italic')
        axis.set_zticks(z_ticks)
        axis.set_zticklabels([0.997, '', 0.998, '', 0.999, '', 1])
        axis.set_zlim(0.997, 1)
        plt.savefig('003_comp_fact(T,P).pdf')
    elif ver == 4: # rho_m
        z_ticks = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
        fig.colorbar(surf, shrink=0.5, aspect=10, pad=0.07, format='%.1f', ticks=z_ticks)
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=0.17, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'$\rho_m$', style='italic')
        axis.set_zticks(z_ticks)
        axis.set_zlim(0.17, 1.46)
        axis.text(322, 19500, 1.35, 'RH = '+str(int(r_hum*100))+'%')
        plt.savefig('004_rho_m(T,P,RH)'+str(int(r_hum*100))+'.pdf')
    elif ver == 5: # m_v
        z_ticks = [0, 0.01, 0.02, 0.03, 0.04, 0.05]
        fig.colorbar(surf, shrink=0.5, aspect=10, pad=0.07, format='%.2f', ticks=z_ticks)
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=-0.008, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'm$_v$', style='italic')
        axis.set_zticks(z_ticks)
        axis.set_zlim(-0.008, 0.055)
        axis.text(285, 100000, 0.04, 'RH = '+str(int(r_hum*100))+'%')
        plt.savefig('005_m_v(T,P,RH)'+str(int(r_hum*100))+'.pdf')
    elif ver == 6: # rho_v
        z_ticks = [0, 0.01, 0.02, 0.03, 0.04, 0.05]
        fig.colorbar(surf, shrink=0.5, aspect=10, pad=0.07, format='%.2f', ticks=z_ticks)
        axis.contour(x_grid, y_grid, z_grid, zdir='z', offset=-0.008, cmap=cm.coolwarm)
        axis.set_zlabel('\n'+r'$\rho_v$, [kg m$^{-3}$]')
        axis.set_zticks(z_ticks)
        axis.set_zlim(-0.008, 0.055)
        axis.text(285, 100000, 0.04, 'RH = '+str(int(r_hum*100))+'%')
        plt.savefig('006_rho_v(T,P,RH)'+str(int(r_hum*100))+'.pdf')


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

def comp_fact(temp_k, p_pa):
    """Use temp_k and p_pa to calculate the compressibility facor

    This function uses the virial equation (14) from Tsilingiris (2008), "Thermophysical and
    transport properties of humid air at temperature range between 0 and 100 C," to calculate the
    compressibility vactor of moist air.

    Positional argument(s):
    temp_k -- temperature in Kelvin
    p_pa   -- total pressure in Pa
    
    Output(s):
    z_c    -- compressibility factor of moist air
    """
    c_coeff = [0.7e-8, -0.147184e-8, 1734.29]
    k_coeff = [0.104e-14, -0.335297e-17, 3645.09]

    a_coeff = c_coeff[0] + c_coeff[1]*exp(c_coeff[2]/temp_k)
    b_coeff = k_coeff[0] + k_coeff[1]*exp(k_coeff[2]/temp_k)

    p_sat = eff_p_sat_liq(temp_k, p_pa)

    return 1 + a_coeff*p_sat + b_coeff*pow(p_sat, 2)

def rho_mix(temp_k, p_pa, r_hum=0.5):
    """Use temp_k, p_pa, and r_hum to calculate the density of moist air

    This function uses the equation of state (1) from  Picard et al. (2008), "Revised formula for
    the density of moist air (CIPM-2007)," to calculate the density of moist air in kg/m^3.

    Positional argument(s):
    temp_k -- temperature in Kelvin
    p_pa   -- total pressure in Pa

    Keyword argument(s):
    r_hum  -- relative humidity as a fraction

    Output(s):
    rho_m  -- density of the mixture in kg/m^3
    """
    z_c = comp_fact(temp_k, p_pa)
    x_v = r_hum*eff_p_sat_liq(temp_k, p_pa)/p_pa
    return (p_pa*28.96546)/(z_c*8314.4598*temp_k)*(1-x_v*(0.3780427))

def x_v(temp_k, p_pa, r_hum=0.5):
    """
    Use temp_k, p_pa, and r_hum to calculate mole fraction

    Positional argument(s):
    temp_k -- temperature in Kelvin
    p_pa   -- total pressure in Pa

    Keyword argument(s):
    r_hum  -- relative humidity as a fraction

    Output(s):
    mole_frac -- mole fraction of the vapor
    """
    return r_hum*eff_p_sat_liq(temp_k, p_pa)/p_pa

def m_v(temp_k, p_pa, r_hum=0.5):
    """
    Use temp_k, p_pa, and r_hum to calculate mass fraction

    Positional argument(s):
    temp_k -- temperature in Kelvin
    p_pa   -- total pressure in Pa

    Keyword argument(s):
    r_hum  -- relative humidity as a fraction

    Output(s):
    mass_frac -- mass fraction of the vapor
    """
    mole_frac = x_v(temp_k, p_pa, r_hum)
    
    return  mole_frac*18.01528/(mole_frac*18.01528 + (1-mole_frac)*28.96546)

def rho_v(temp_k, p_pa, r_hum=0.5):
    """
    Use temp_k, p_pa, and r_hum to calculate density of the vapor

    Positional argument(s):
    temp_k -- temperature in Kelvin
    p_pa   -- total pressure in Pa

    Keyword argument(s):
    r_hum  -- relative humidity as a fraction

    Output(s):
    rho_v -- density of the vapor
    """
    return m_v(temp_k, p_pa, r_hum)*rho_mix(temp_k, p_pa, r_hum)

def main():
    """Main"""
