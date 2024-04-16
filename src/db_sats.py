import numpy as np
import ufl

def S_HbXY(P_O2, P_CO2, dash_params):
    '''
    Saturations S_HbO2 and S_HbCO2 as functions of
    gas partial pressures P_O2 and P_CO2. It is assumed
    that partial pressures are equal in plasma and inside
    red blood cells.

    Functions defined in Dash, Bassingthwaighte & Korman articles:
    - R. K. Dash, B. Korman, J. B. Bassingthwaighte. (2016). Simple accurate mathematical models of blood HbO2 and HbCO2 dissociation curves at
    varied physiological conditions: Evaluation and comparison with other models, European Journal of Applied Physiology 116 (1)
    97-113. doi:10.1007/s00421-015-3228-3.
    - R. K. Dash, J. B. Bassingthwaighte. (2010). Erratum to: Blood HbO2 and HbCO2 Dissociation Curves at Varied O2, CO2, pH, 2,3-DPG and
    Temperature Levels, Annals of Biomedical Engineering 38 (4) 1683-1701. doi:10.1007/s10439-010-9948-y.
    '''

    # Known parameters (reaction rates, Fick's law constants)
    alpha_CO2 = dash_params["alpha_CO2"]        # M mmHg^-1
    alpha_O2 = dash_params["alpha_O2"]          # M mmHg^-1

    K2p = dash_params["K2p"]                    # M^-1
    K2pp = dash_params["K2pp"]                  # M
    K3p = dash_params["K3p"]                    # M^-1
    K3pp = dash_params["K3pp"]                  # M
    K5pp = dash_params["K5pp"]                  # M
    K6pp = dash_params["K6pp"]                  # M

    # Red blood cell pH (with pH_plasma = 7.4).
    pH_rbc = dash_params["pH_rbc"] 
    pH_pla = dash_params["pH_pla"]

    R_rbc = 10**(-pH_pla+pH_rbc)
    c_H = 10**(-pH_rbc)

    # Variables defined in Dash & Bassingthwaighte (2016).
    phi_1 = 1 + (K2pp/c_H)
    phi_2 = 1 + (K3pp/c_H)
    phi_3 = 1 + (c_H/K5pp)
    phi_4 = 1 + (c_H/K6pp)

    # P50 as function of carbon dioxide deviation from standard.
    P50_std = dash_params["P50_std"] 
    P_CO2_std = dash_params["P_CO2_std"] 
    pH_rbc_S = dash_params["pH_rbc_S"]
    P50_delta_CO2 = P50_std + 1.273E-1*(P_CO2-P_CO2_std) + 1.083E-4*(P_CO2-P_CO2_std)**2
    P50_delta_pH = P50_std - 25.535*(pH_rbc - pH_rbc_S) + 10.646*(pH_rbc - pH_rbc_S)**2 -1.764*(pH_rbc - pH_rbc_S)**3
    P50 = (P50_delta_CO2*P50_delta_pH)/P50_std

    # Variable nH (Hill equation exponent), fit to data.
    alpha = dash_params["alpha"]
    beta = dash_params["beta"]
    gamma = dash_params["gamma"]
    nH = alpha - beta*10**(-P_O2/gamma)

    # Constant K4'
    K4p = ((alpha_O2*P_O2)**(nH-1))*(K2p*alpha_CO2*P_CO2*phi_1 + phi_3)
    K4p /= (((alpha_O2*P50)**(nH))*(K3p*alpha_CO2*P_CO2*phi_2 + phi_4))

    # Variable parameters for both saturations
    K_HbO2 = (K4p*(K3p*phi_2*alpha_CO2*P_CO2 + phi_4))
    K_HbO2 /= (K2p*phi_1*alpha_CO2*P_CO2 + phi_3)

    K_HbCO2 = (K2p*phi_1 + K3p*K4p*phi_2*alpha_O2*P_O2)
    K_HbCO2 /= (phi_3 + K4p*phi_4*alpha_O2*P_O2)

    # Final saturations
    S_HbO2 = ufl.variable((K_HbO2*alpha_O2*P_O2)/(1+(K_HbO2*alpha_O2*P_O2)))
    S_HbCO2 = ufl.variable((K_HbCO2*alpha_CO2*P_CO2)/(1+(K_HbCO2*alpha_CO2*P_CO2)))
    
    return S_HbO2, S_HbCO2