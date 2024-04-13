def instance_params():
    perfusion_params = {
        "mu": 22.5E-6,                    # mmHg s
        "kappa": 0.13126,                 # um^2
        "pmin": 8.0,                      # mmHg
        "uin": 160.0                      # um/s
    }

    dash_params = {
        "alpha_CO2": 3.27E-5,             # M mmHg^-1
        "alpha_O2": 1.46E-6 ,             # M mmHg^-1
        "K2p": 21.5,                      # M^-1
        "K2pp": 1E-6,                     # M
        "K3p": 11.3,                      # M^-1
        "K3pp": 1E-6,                     # M
        "K5pp": 2.63E-8,                  # M
        "K6pp": 1.91E-8,                  # M
        "pH_rbc": 7.24,                   # Red blood cell pH (with pH_plasma = 7.4).
        "pH_rbc_S": 7.24,
        "pH_pla": 7.4,                    # plasma pH
        "P50_std": 26.8,                  # mmHg
        "P_CO2_std": 40,                  # mmHg
        "alpha": 2.82,                    # parameters for nH calculation
        "beta": 1.20,
        "gamma": 29.25,
        "W_bl": 0.81,                     # nondimensional, con Hct = 0.45
        "beta_O2": 1.46E-6*1E-12,         # mmol/(um3*mmHg)
        "beta_CO2": 3.27E-5*1E-12,        # mmol/(um3*mmHg)
        "Hb_bl": 2.33E-3*1E-12,           # mmol/um3
        "R_rbc": 0.692,                   # nondimensional, con pH_pl = 7.4
        "Hct": 0.45,                      # nondimensional (volumetric rbc fraction in blood)
        "W_rbc": 0.65,                    # nondimensional
        "W_pl": 0.94,                     # nondimensional
        "K1": 7.94E-7*1E-12               # mmol/um3 (hydration constant CO2 + H2O = HCO3 + H)
    }

    transport_params = {
        "d_pla_O2": 1.62E3,               # um2/s
        "d_pla_CO2": 1E3,                 # um2/s
        "d_ba_O2": 0.914E3,               # um2/s
        "d_ba_CO2": 1E3,                  # um2/s
        "h_ba": 0.3,                      # um
        "p_O2_air": 100,                  # mmHg
        "p_CO2_air": 40,                  # mmHg
        "p_O2_in": 40,                    # mmHg
        "p_CO2_in": 45                    # mmHg
    }

    return perfusion_params, dash_params, transport_params