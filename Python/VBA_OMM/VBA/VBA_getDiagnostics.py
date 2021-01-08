import numpy as np
from . import VBA_basics as base


def get_fit(data, priors, posterior, options, suffStat):

    dim = options["dim"]

    # No. of free Parameters
    n_p = dim["n_phi"] + dim["n_theta"] + dim["n"]

    # Log-Likelihood
    v = posterior["b"]/posterior["a"]
    LL = -0.5 * suffStat["dy2"] / v
    for i in range(0, dim["nD"]):
        ldQ = base.log_det(priors["iQy"][i]/v)
        LL = LL + 0.5*ldQ
    LL = LL - 0.5 * dim["nD"] * np.log(2*np.pi)

    # AIC/BIC
    AIC = LL - n_p
    BIC = LL - 0.5 * n_p * np.log(dim["nD"])

    # R2
    y = np.reshape(data["y"], dim["nD"])
    R2 = 1 - (suffStat["dy2"] / np.sum((y - np.mean(y)) ** 2))
    if R2<0:
        R2 = 0

    fit = {"np": n_p,
           "LL": float(LL),
           "BIC": float(BIC),
           "AIC": float(AIC),
           "R2": float(R2)}

    return fit
