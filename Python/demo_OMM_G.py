import numpy as np
import VBA_OMM
import pandas as pd

# Read the demo data from csv file
df = pd.read_csv("DataOMM_G.csv")
t = df["Time [min]"].to_numpy()
G = df["Glucose [mmol/L]"].to_numpy()
I = df["Insulin [pmol/L]"].to_numpy()

# Construt the data struture
dat = {"t": t.reshape(len(t),1),
       "G": G.reshape(len(t),1),
       "I": I.reshape(len(t),1)}

# Constants
const = {"A": 6,                    # ??
         "V": 0.145,                # L/kg
         "dt": 0.1,                 # min
         "Rap": [],                 # mmol/kg/min
         "X0": 0,                   # 1/min
         "measCV": 2,               # %
         "Gb": dat["G"][0, 0],      # mmol/L
         "G0": dat["G"][0, 0],      # mmol/L
         "Ib": dat["I"][0, 0]}      # pmol/L

# Construct inversion options
opt = {"GA_fun": "RaPL",
       "tb": np.array([0, 10, 30, 60, 90, 120, 180, 300]),      # min
       "alpha": 0.017,                                          # 1/min
       "displayWin": True}

# Priors
    # - System Parameters [median CV]
priors = {"p1": np.array([0.025, 25]),      # 1/min
          "p2": np.array([0.012, 40]),      # 1/min
          "SI": np.array([12E-5, 100])}     # 1/min per pmol/L

    # - Input function Parameters
if opt["GA_fun"] == 'RaPL':
    priors.update({"k": np.array([[3.2E-3*const["A"], 50],          # mmol/kg/min
                                  [7.3E-3*const["A"], 50],          # mmol/kg/min
                                  [5.4E-3*const["A"], 50],          # mmol/kg/min
                                  [5.1E-3*const["A"], 50],          # mmol/kg/min
                                  [3.7E-3*const["A"], 50],          # mmol/kg/min    
                                  [1.8E-3*const["A"], 50]])})       # mmol/kg/min
if opt["GA_fun"] == 'RaLN':
    priors.update({"k": np.array([[30, 30],         # min
                                  [0.5, 30],        
                                  [100, 30],        # min
                                  [0.5, 30],
                                  [0.7, 30]])})

out = VBA_OMM.mainG(dat, priors, const, opt)
