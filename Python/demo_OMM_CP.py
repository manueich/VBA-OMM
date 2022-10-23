import numpy as np
import VBA_OMM
import pandas as pd

# Read the demo data from csv file
df = pd.read_csv("Data2.csv")
t = df["Time [min]"].to_numpy()
G = df["Glucose [mmol/L]"].to_numpy()
CP = df["C-Peptide [nmol/L]"].to_numpy()

# Construt the data struture
dat = {"t": t.reshape(len(t),1),
       "G": G.reshape(len(t),1),
       "CP": CP.reshape(len(t),1)}

# Constants
const = {"dt": 0.1,
         "measCV": 6,
         "age": 37,
         "subject_type": "obese",
         "CPb": dat["CP"][0, 0],
         "Gb": dat["G"][0, 0]}

# Construct inversion options
opt = {"displayWin": True,
       "updateMeasCV": False}

# Priors [median CV]
priors = {"T": np.array([10, 50]),                 # min 
          "beta": np.array([20, 50]),              # 1E-9 1/min
          "h": np.array([dat["G"][0,0], 20]),      # mmol/l
          "kd": np.array([1000, 50])}              # 1E-9

out = VBA_OMM.mainCP(dat, priors, const, opt)
