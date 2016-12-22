"""
log and exp used for evaluating formulas
"""
from math import log, exp

"""This python modual containts functions for the calculation of thermophysical properties of moist
air in the range between 0 and 100 degC. For details regarding these functions see P.T. Tsilingiris
(2008) "Thermophsical and transport properties of humid air at temperature range between 0 and 100
degC," Energy Conv. Manag. (49) pp. 1098-1110.

Author: Rich H. Inman
email: rinman@ucsd.edu
Date Created: 21 Nov. 2016
"""

def check_temp_liq(temp_k, delta=0):
    """Check the temperature in Kelvin to return boolean of liquid water

    This funciton checkes the temperature in Kelvin to see if the temperature is between 273.15 and
    373.15 K. Returns True if the temperature is within said range and False otherwise.

    Input(s):
        temp_k -> temperature in Kelvin
        delta  -> optional uncertainty in the measurement of temp_k (default is zero).
                      This input is used when calculating uncertainties and is added to the logical
                      checks of the temperature range. For example, if you want to calculate the
                      uncertainty at 373.15 K, which is the upper bound for this function, you would
                      need to evaluate this function at 373.15 K + delta.
    Output(s):
        bool   -> True if the temperature is between 273.15 and 373.15 K else False
    """
    if temp_k < 273.15:
        print '\tError in p_sat_l: Temperature below 273.15 K.'
        return False
    elif temp_k > 373.15 + delta:
        print '\tError in p_sat_l: Temperature above 373.15 K.'
        return False
    else:
        return True

def p_sat_liq(temp_k):
    """Use Hardy (1998) to return P_sat(T_K) in Pa.

    This function uses Hardy (1998), "ITS-90 Formulations for Vapor Pressure, Frostpoint
    Temperature, Dewpoint Temperature, and Enhancement Factors in the Range -100 to +100 C," to
    calculate the saturation vapor pressure over liquid water in the range 0 to 100 C.

    Input(s):
        temp_k -> temperature in Kelvin
    Output(s):
        p_sat  -> saturation vapor pressure over liquid water in Pa
                      False if the input temp_k is outside of the limits 272.15 to 353.15 + delta
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

    Input(s):
        temp_k -> temperature in Kelvin
        p_pa   -> total pressure in Pa
        p_sat  -> saturation vapor pressure of water over liquid water in Pa
    Output(s):
        f_enh -> enhancement factor over liquid water
                      False if the input temp_k is outside of the limits 272.15 to 353.15 + delta
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

    Input(s):
        temp_k -> temperature in Kelvin
        p_pa   -> total pressure in Pa
    Output(s):
        eff_p  -> effective saturation vapor pressure over liquid water
                      False if the input temp_k is outside of the limits 272.15 to 353.15 + delta
    """
    p_sat = p_sat_liq(temp_k)
    enh_fact = enh_fact_liq(temp_k, p_pa, p_sat)
    return enh_fact*p_sat

def get_eff_p_sat_liq_delta(temp_k, p_pa, delta=0.2):
    """Use eff_p_sat_liq twice to return eff_p_sat with uncertainty.

    This function uses  Hardy (1998), "ITS-90 Formulations for Vapor Pressure, Frostpoint
    Temperature, Dewpoint Temperature, and Enhancement Factors in the Range -100 to +100 C," to
    calculate the effective saturation vapor pressure over liquid water in the range 0 to 100 C as
    well as the associated uncertainty. The uncertainty in effective saturation vapor pressure is
    propigated from the uncertainty in the measurement of the temperature. The default uncertainty
    is 0.2 K. If you would like to chage it, use the optional delta parameter; e.g., delta=0.1. It
    should be noted that there is also an uncertainty associated with the pressure measurement, but
    because the saturation vapor pressure is a strong function of temperature and a only a weak
    function of pressure, the uncertainty from the pressure is negligible.

    Input(s):
        temp_k -> temperature in Kelvin
        p_pa   -> total pressure in Pa
    Optional Input(s):
        delta  -> uncertainty in the measurement of the temperature
    Output(s):
        (p_sat, delta_p_sat) -> a tuple of the saturation vapor pressure and its uncertainty in Pa.
    """
    p_1 = eff_p_sat_liq(temp_k, p_pa)
    p_2 = eff_p_sat_liq(temp_k + delta, p_pa)
    return (p_1, abs(p_2 - p_1))

def main():
    """Main"""
