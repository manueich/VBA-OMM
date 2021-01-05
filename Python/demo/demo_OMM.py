import numpy as np
from Functions import VBA_OMM_G as Inv
import csv

# Read the data from csv file
with open('demo_dat.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    t = np.zeros((1, 19))
    G = np.zeros((1, 19))
    I = np.zeros((1, 19))
    i = 0
    for row in csv_reader:
        t[0, i] = row[0]
        G[0, i] = row[1]
        I[0, i] = row[2]
        i=i+1

# Construt the data struture
dat = {"t": t,
       "G": G,
       "I": I}

# Constants
const = {"A": 6,
         "V": 0.145,
         "dt": 1,
         "Rap": [],
         "X0": 0,
         "measCV": 2,
         "Gb": dat["G"][0, 0],
         "G0": dat["G"][0, 0],
         "Ib": dat["I"][0, 0]}

# Construct inversion options
opt = {"GA_fun": "RaLN",
       "tb": np.array([0, 10, 30, 60, 90, 120, 180, 300]),
       "alpha": 0.017,
       "displayWin": True}

# Priors
    # - System Parameters [median CV]
priors = {"p1": np.array([0.025, 25]),
          "p2": np.array([0.012, 40]),
          "SI": np.array([7.1E-4, 100])}

    # - Input function Parameters
if opt["GA_fun"] == 'RaPL':
    priors.update({"k": np.array([[3.2E-3*const["A"], 50],
                                  [7.3E-3*const["A"], 50],
                                  [5.4E-3*const["A"], 50],
                                  [5.1E-3*const["A"], 50],
                                  [3.7E-3*const["A"], 50],
                                  [1.8E-3*const["A"], 50]])})
if opt["GA_fun"] == 'RaLN':
    priors.update({"k": np.array([[30, 30],
                                  [0.5, 30],
                                  [100, 30],
                                  [0.5, 30],
                                  [0.7, 30]])})

out = Inv.VBA_OMM_G(dat, priors, const, opt)

print("ja")

