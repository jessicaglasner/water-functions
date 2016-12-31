"""
log and exp used for evaluating formulas
errorAnalysis for error analysis
"""
from math import log, exp
import errorAnalysis as err

"""This python modual containts functions for the calculation of thermophysical properties of moist
air in the range between 0 and 100 degC. For details regarding these functions see P.T. Tsilingiris
(2008) "Thermophsical and transport properties of humid air at temperature range between 0 and 100
degC," Energy Conv. Manag. (49) pp. 1098-1110.

Author: Rich H. Inman
email: rinman@ucsd.edu
Date Created: 21 Nov. 2016
"""

def get_temp_k_obs(temp_k, delta_t=0.2):
    """Use temp k to return list of tuples containing uncertainties.

    This function accepts either a single observation or an iterable of observations and appends the
    uncertainty to each of the observations. The default uncertainty is 0.2 K.

    Positional argument(s):
    temp_k  -- temperature(s) in Kelvin

    Keyword argument(s):
    delta_t -- uncertainty of all temperature observations (default 0.2)

    Output(s)
    temp_k_obs -- list of tuples [(obs_1, delta_1), ..., (obs_n, delta_n)]
    """
    try:
        temp_k_obs = [(temp, delta_t) for temp in temp_k]
    except TypeError:
        temp_k_obs = [(temp_k, delta_t)]
    return temp_k_obs

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
