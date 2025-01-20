import numpy as np
import scipy.linalg as la
import scipy.special as spec
import warnings


def solveODE(t, muP, SigmaP, u,  options):
    # Function solves the ODE and provides the sensitivity equations,
    # i.e. derivatives of model states and output with respect to parameters and initial conditions

    dim = options["dim"]

    phi = muP[0: dim["n_phi"]]
    th = muP[dim["n_phi"]: dim["n_theta"]+dim["n_phi"]]
    x0 = muP[dim["n_theta"]+dim["n_phi"]: dim["n_theta"]+dim["n_phi"]+dim["n"]]

    # Get functions
    f_model = options["f_model"]
    f_obs = options["f_obs"]

    # Preallocate Variables
    x = np.zeros((dim["n"], np.shape(t)[1]))
    y = np.zeros((dim["nY"], np.shape(t)[1]))
    SigmaX = [SigmaP[dim["n_theta"]+dim["n_phi"]: dim["n_theta"]+dim["n_phi"]+dim["n"], dim["n_theta"]+dim["n_phi"]: dim["n_theta"]+dim["n_phi"]+dim["n"]]]

    # Set initial conditions
    for i in range(0, dim["n"]):
        x[i, 0] = x0[i]
    y[:, [0]] = f_obs(x[:, 0], phi, u[:, 0], options["inG"])[0]

    # Preallocate Variables
    dXdTh = [np.zeros((dim["n"], dim["n_theta"]))]
    for i in range(0, np.shape(t)[1] - 1):
        dXdTh.append(np.zeros((dim["n"], dim["n_theta"])))

    dXdX0 = [np.eye(dim["n"])]
    for i in range(0, np.shape(t)[1] - 1):
        dXdX0.append(np.zeros((dim["n"], dim["n"])))

    dYdTh = [f_obs(x[:, 0], phi, u[:, 0], options["inG"])[1] @ dXdTh[0]]
    dYdX0 = [f_obs(x[:, 0], phi, u[:, 0], options["inG"])[1] @ dXdX0[0]]
    dYdPhi = [f_obs(x[:, 0], phi, u[:, 0], options["inG"])[2]]

    dG_dP = [np.concatenate((dYdPhi[0].T, dYdTh[0].T, dYdX0[0].T), 0)]

    dt = t[0, 1] - t[0, 0]

    # Loop over integration time
    for i in range(0, np.shape(t)[1] - 1):

        # Model
        if options["ODESolver"] == 'RK':
            x[:, [i + 1]], dXdTh[i + 1], dXdX0[i + 1] = Runge_Kutta(f_model, i, x[:, [i]], dXdTh[i], dXdX0[i], u, th, options, dt)
        if options["ODESolver"] == 'Euler':
            x[:, [i + 1]], dXdTh[i + 1], dXdX0[i + 1] = Euler(f_model, i, x[:, [i]], dXdTh[i], dXdX0[i], u, th, options, dt)
        # Observation
        y[:, [i + 1]], dY_dX, dY_dPhi = f_obs(x[:, [i + 1]], phi, u[:, [i]], options["inG"])

        # Calculate Sensitivities of output wrt
        dYdTh.append(dY_dX @ dXdTh[i + 1])   # Evolution Parameters
        dYdX0.append(dY_dX @ dXdX0[i + 1])   # Initial Conditions
        dYdPhi.append(dY_dPhi)               # Observation Parameters
        dG_dP.append(np.concatenate((dYdPhi[i + 1].T, dYdTh[i + 1].T, dYdX0[i + 1].T), 0))

        # Uncertainties in model states
        dXdP = np.concatenate((dXdTh[i+1], dXdX0[i+1]), 1)
        SigmaX.append(dXdP @ SigmaP[dim["n_phi"]: dim["n_phi"]+dim["n_theta"]+dim["n"], dim["n_phi"]: dim["n_phi"]+dim["n_theta"]+dim["n"]] @ dXdP.T)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

    return y, x, SigmaX, dXdTh, dXdX0, dY_dPhi, dYdTh, dYdX0, dG_dP


def Euler(f_model, i, xn, Sthn, Sx0n, u, th, options, dt):

    un = u[:, [i]]

    dx, dFdX, dFdth = f_model(xn, th, un, options["inF"])
    xn2 = xn + dt * dx
    Sthn2 = Sthn + dt * (dFdth + dFdX @ Sthn)
    Sx0n2 = Sx0n + dt * (dFdX @ Sx0n)

    return xn2, Sthn2, Sx0n2


def Runge_Kutta(f_model, i, xn, Sthn, Sx0n, u, th, options, dt):

    un = u[:, [i]]
    un2 = u[:, [i+1]]

    # Model States
    k1 = f_model(xn, th, un, options["inF"])[0]
    k2 = f_model(xn + dt*k1/2, th, (un + un2)/2, options["inF"])[0]
    k3 = f_model(xn + dt*k2/2, th, (un + un2)/2, options["inF"])[0]
    k4 = f_model(xn + dt*k3, th, un2, options["inF"])[0]

    xn2 = xn + dt/6*(k1 + 2*k2 + 2*k3 + k4)

    # Calculate Sensitivities of model states wrt
        # Evolution Parameters

    dx, dFdX, dFdth = f_model(xn, th, un, options["inF"])
    k1 = dFdth + dFdX @ Sthn
    k2 = dFdth + dFdX @ (Sthn + dt * k1 / 2)
    k3 = dFdth + dFdX @ (Sthn + dt * k2 / 2)
    k4 = dFdth + dFdX @ (Sthn + dt * k3)

    Sthn2 = Sthn + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

        # Initial Conditions
    k1 = dFdX @ Sx0n
    k2 = dFdX @ (Sx0n + dt * k1 / 2)
    k3 = dFdX @ (Sx0n + dt * k2 / 2)
    k4 = dFdX @ (Sx0n + dt * k3)

    Sx0n2 = Sx0n + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    return xn2, Sthn2, Sx0n2


def Free_Energy(posterior, priors, suffStat, options):
    # Calculates the Free Energy

    # Entropy Calculus
        # System Parameters
    SP = 0.5 * (options["dim"]["n"] + options["dim"]["n_theta"] + options["dim"]["n_phi"]) * np.log(2 * np.pi * np.exp(1)) + 0.5 * log_det(posterior["SigmaP"])
     #   # Noise Parameters
    # Ssigma = spec.gammaln(posterior["a"]) - np.log(posterior["b"]) + (1 - posterior["a"]) * spec.digamma(
        # posterior["a"]) + posterior["a"]

    # Noise
    SSE = 0
    dF = 0
    ldQ = 0
    E = posterior["a"]/posterior["b"]
    V = posterior["a"]/(posterior["b"]**2)
    E0 = priors["a"]/priors["b"]
    V0 = priors["a"]/(priors["b"]**2)
    dF = dF - D_KL(E, V, E0, V0, 'Gamma')
    SSE = SSE + E * suffStat["dy2"]
    ElogS = spec.digamma(posterior["a"]) - np.log(posterior["b"])

    dF = dF + 0.5*options["dim"]["nD"]*options["dim"]["nY"]*ElogS
    for i in range(0, options["dim"]["nD"]):
        ldQ = ldQ + log_det(priors["iQy"][i])

    # Parameters
    ntot = (options["dim"]["nD"] * options["dim"]["nY"]) + options["dim"]["n"] + options["dim"]["n_theta"] + options["dim"]["n_phi"]
    SSE = SSE + suffStat["dP"].T @ Invert_M(priors["SigmaP"]) @ suffStat["dP"]
    ldQ = ldQ - log_det(priors["SigmaP"])
    S = SP - 0.5*(options["dim"]["n"] + options["dim"]["n_theta"] + options["dim"]["n_phi"])

    # Compose Free Energy
    F = -0.5*SSE - 0.5*ntot*np.log(2*np.pi) + 0.5*ldQ + S + dF

    return F


def log_det(Q):
    # Computes the log determinant of matrix Q

    tol = np.spacing(1)     # Floating-Point Tolerance

    # Check Matrix
    if np.allclose(Q, np.eye(Q.shape[0]), rtol=1e-12, atol=1e-12):  # Identity Matrix
        ldQ = 0
    elif np.allclose(Q, np.diag(np.diag(Q)), rtol=1e-12, atol=1e-12):  # Diagonal Matrix
        dQ = np.diag(Q)
        ldQ = 0
        for i in range(0, dQ.size):
            if np.abs(dQ[i]) > tol:
                ldQ = ldQ + np.log(dQ[i])
    else:  # Full Matrix
        if not isWeird(Q):
            try:
                ldQ = np.log(la.det(Q))
            except:
                try:
                    ldQ = np.linalg.slogdet(Q)[1]
                except:
                    raise Exception('ERROR: Calculation of log determinant of matrix has failed')
        else:
            ldQ = 0

    return np.real(ldQ)


def Invert_M(Q):
    # Invert Matrix possibly containing NaNs and Infs

    dq = np.diag(Q)
    # Detect NaNs or Infs on the main diagonal
    idx = np.where(np.logical_and(np.logical_and(np.invert(np.isnan(dq)), dq != 0), np.invert(np.isinf(dq))))

    if np.size(idx) == 0:   # All NaNs -> return zero matrix
        iQ = np.zeros((Q.shape[0], Q.shape[0]))
        return iQ
    elif np.size(idx) == Q.shape[0]:    # No NaNs
        issub = False
        subQ = Q
    else:       # Some NaNs -> extract sub matrix from non-NaN entries
        issub = True
        k = 0
        subQ = np.zeros((np.size(idx)**2))
        for i in range(0, Q.shape[0]):
            for j in range(0, Q.shape[0]):
                if np.any(np.isin(idx, i)) and np.any(np.isin(idx, j)):
                    subQ[k] = Q[i, j]
                    k = k+1
        subQ = subQ.reshape(np.size(idx), np.size(idx))

    # Invert matrix
    if np.allclose(subQ, np.eye(subQ.shape[0]), rtol=1e-12, atol=1e-12):    # Identity Matrix
        subiQ = subQ
    elif np.allclose(subQ, np.diag(np.diag(subQ)), rtol=1e-12, atol=1e-12):    # Diagonal Matrix
        # Find Tolerance of float precision
        nrm = np.spacing(la.norm(np.diag(subQ), np.inf)) * subQ.shape[0]
        if nrm <= np.exp(-32):
            tol = np.exp(-32)
        else:
            tol = nrm
        subiQ = np.diag(1/(np.diag(subQ)+tol))
    else:                                                                      # Full Matrix
        # Find Tolerance of float precision
        nrm = np.spacing(la.norm(np.diag(subQ), np.inf)) * subQ.shape[0]
        if nrm <= np.exp(-32):
            tol = np.exp(-32)
        else:
            tol = nrm
        subiQ = la.inv(subQ + np.eye(subQ.shape[0]) * tol)

    # Reassemble to full sized matrix if NaNs were contained
    if issub:
        iQ = np.zeros((Q.shape[0], Q.shape[0]))
        subiQ = subiQ.reshape((1, np.size(idx)**2))
        k = 0
        for i in range(0, iQ.shape[0]):
            for j in range(0, iQ.shape[0]):
                if np.any(np.isin(idx, i)) and np.any(np.isin(idx, j)):
                    iQ[i, j] = subiQ[0, k]
                    k = k+1
    else:
        iQ = subiQ

    return iQ


def D_KL(m1, v1, m2, v2, dist):
    # Computes the Kullback-Leibler divergence between distributions (Normal or Gamma) with
    # first and second order moments m and v

    if np.all(m1 == m2) and np.all(v1 == v2):
        DKL = 0
        return DKL

    if dist == "Normal":
        n = np.size(m1, 0)
        iv2 = Invert_M(v2)
        vv = v1 @ iv2
        DKL = 0.5*(-log_det(vv) + np.trace(vv) + (m1-m2).T @ iv2 @ (m1-m2) - n)
        DKL = float(DKL)

    if dist == "Gamma":
        b1 = m1/v1
        b2 = m2/v2
        a1 = b1*m1
        a2 = b2*m2

        DKL = spec.gammaln(a2) - spec.gammaln(a1) - a1*(1-b2/b1) + a2*np.log(b1/b2) + (a1-a2)*spec.digamma(a1)

    return DKL


def isWeird(X):
    # Check if input array contains NaNs, Infs or complex numbers

    if np.any(np.isnan(X)) or np.any(np.isinf(X)) or np.any(np.iscomplex(X)):
        return True
    else:
        return False

