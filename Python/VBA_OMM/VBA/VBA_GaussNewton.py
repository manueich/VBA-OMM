import numpy as np
from . import VBA_basics as base
import warnings


def GaussNewton(data, t, posterior, priors, suffStat, options):
    # Initialise
    dim = options["dim"]

    previousMu = np.concatenate((posterior["muPhi"], posterior["muTheta"], posterior["muX0"]), 0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        previousI, dmu, Sigma, suffStat2 = get_move(previousMu, data, t, posterior, priors, suffStat.copy(), options)

    # Gauss-Newton Update routine
    stop = False
    it = 0
    conv = False
    posterior0 = posterior.copy()
    posterior.update({"muP": previousMu})
    posterior.update({"SigmaP": Sigma})

    while not stop:
        it = it + 1

        # Make move
        mu = previousMu + dmu

        # Get next move
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            I, NextdeltaMu, Sigma, suffStat2 = get_move(mu, data, t, posterior, priors, suffStat.copy(), options)
        deltaI = I - previousI
        if previousI == 0:
            rdf = np.NaN
        else:
            rdf = deltaI / np.abs(previousI)

        # Accept move or halve step size
        if deltaI < 0:  # Half step size
            dmu = dmu / 2
        else:  # Accept mu
            dmu = NextdeltaMu
            previousMu = mu
            previousI = I

            # Update Posterior
            update_posterior(posterior, dim, previousMu, Sigma)

            # Calculate new Free Energy
            F = base.Free_Energy(posterior, priors, suffStat2, options)
            Fall = suffStat2["F"]
            Fall.append(F)

            # Update sufficient statistics
            suffStat2.update({"F": Fall})
            suffStat = suffStat2.copy()

            conv = True

        # Check convergence criterion
        if np.abs(rdf) <= options["GnTolFun"] or it == options["GnMaxIter"]:
            stop = True

    if not conv:
        Fall = suffStat["F"]
        Fall.append(Fall[-1])
        suffStat.update({"F": Fall})
        posterior = posterior0.copy()

    return posterior, suffStat


def get_move(mu, data, t, posterior, priors, suffStat, options):
    sigmaHat = posterior["a"] / posterior["b"]

    yd = data["y"]
    td = data["t"]

    # Preallocate variables
    Q = priors["SigmaP"]
    iQ = base.Invert_M(Q)
    mu0 = priors["muP"]
    dmu0 = mu0 - mu
    ddYdP = 0
    dy = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))
    vy = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))
    gx = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))
    dy2 = 0
    d2GdX2 = 0
    div = False

    # Get estimates at current mode
    y, muX, SigmaX, dXdTh, dXdX0, dYdPhi, dYdTh, dYdX0, dG_dP = base.solveODE(t, mu, posterior["SigmaP"], data["u"], options)

    # Loop over time series
    for i in range(0, td.size):
        idx = np.where(np.round(t, 8) == np.round(td[0, i], 8))
        gx[:, [i]] = y[:, idx[1]]
        iQyt = priors["iQy"][i]

        # Posterior Covaraince Matrix terms
        d2GdX2 = d2GdX2 + dG_dP[int(idx[1])] @ iQyt @ dG_dP[int(idx[1])].T

        # Error terms
        dy[:, [i]] = yd[:, [i]] - gx[:, [i]]
        dy2 = dy2 + dy[:, [i]].T @ iQyt @ dy[:, [i]]
        ddYdP = ddYdP + dG_dP[int(idx[1])] @ iQyt @ dy[:, [i]]

        # Predictive density (data space)
        V = dG_dP[int(idx[1])].T @ posterior["SigmaP"] @ dG_dP[int(idx[1])] + 1 / sigmaHat * base.Invert_M(iQyt)
        vy[:, [i]] = (np.diag(V)).reshape((options["dim"]["nY"], 1))

        if base.isWeird(dy2) or base.isWeird(dG_dP[int(idx[1])]):
            div = True
            break

    # Posterior Covariance Matrix
    iSigma = iQ + sigmaHat * d2GdX2
    Sigma = base.Invert_M(iSigma)

    # Move
    dmu = Sigma @ (iQ @ dmu0 + sigmaHat * ddYdP)

    # Variational Energy
    I = -0.5 * (dmu0.T @ iQ @ dmu0) - 0.5 * sigmaHat * dy2
    if base.isWeird(I) or base.isWeird(Sigma) or div:
        I = -np.inf

    # Update suffstat
    model_out = {"t":t,
                 "y": y,
                 "muX": muX,
                 "SigmaX": SigmaX,
                 "dXdTh": dXdTh,
                 "dXdX0": dXdX0,
                 "dYdPhi": dYdPhi,
                 "dYdTh": dYdTh,
                 "dYdX0": dYdX0,
                 "dG_dP": dG_dP}

    suffStat.update({"Iphi": I,
                     "gx": gx,
                     "dy": dy,
                     "dy2": dy2,
                     "vy": vy,
                     "dP": dmu0,
                     "div": div,
                     "model_out": model_out,
                     "SigmaP": posterior["SigmaP"]})

    return I, dmu, Sigma, suffStat


def update_posterior(posterior, dim, previousMu, Sigma):
    # Update muP and SigmaP
    posterior.update({"muP": previousMu})
    posterior.update({"SigmaP": Sigma})

    # Disentangle muP and SigmaP
    if dim["n_phi"] == 0 and dim["n_theta"] == 0:
        posterior.update({"muX0": previousMu})
        posterior.update({"SigmaX0": Sigma})
    elif dim["n_phi"] == 0:
        posterior.update({"muTheta": previousMu[0: dim["n_theta"]]})
        posterior.update({"muX0": previousMu[dim["n_theta"]: dim["n_theta"] + dim["n"]]})
        posterior.update({"SigmaTheta": Sigma[0: dim["n_theta"],
                                        0: dim["n_theta"]]})
        posterior.update({"SigmaX0": Sigma[dim["n_theta"]: dim["n_theta"] + dim["n"],
                                     dim["n_theta"]: dim["n_theta"] + dim["n"]]})
    elif dim["n_theta"] == 0:
        posterior.update({"muPhi": previousMu[0: dim["n_phi"]]})
        posterior.update({"muX0": previousMu[dim["n_phi"]: dim["n_phi"] + dim["n"]]})
        posterior.update({"SigmaPhi": Sigma[0: dim["n_phi"],
                                        0: dim["n_phi"]]})
        posterior.update({"SigmaX0": Sigma[dim["n_phi"]: dim["n_phi"] + dim["n"],
                                     dim["n_phi"]: dim["n_phi"] + dim["n"]]})
    else:
        posterior.update({"muPhi": previousMu[0: dim["n_phi"]]})
        posterior.update({"muTheta": previousMu[dim["n_phi"]: dim["n_theta"] + dim["n_phi"]]})
        posterior.update({"muX0": previousMu[dim["n_theta"] + dim["n_phi"]: dim["n_theta"] + dim["n_phi"] + dim["n"]]})
        posterior.update({"SigmaPhi": Sigma[0: dim["n_phi"],
                                      0: dim["n_phi"]]})
        posterior.update({"SigmaTheta": Sigma[dim["n_phi"]: dim["n_theta"] + dim["n_phi"],
                                        dim["n_phi"]: dim["n_theta"] + dim["n_phi"]]})
        posterior.update({"SigmaX0": Sigma[dim["n_theta"] + dim["n_phi"]: dim["n_theta"] + dim["n_phi"] + dim["n"],
                                     dim["n_theta"] + dim["n_phi"]: dim["n_theta"] + dim["n_phi"] + dim["n"]]})

    return posterior
