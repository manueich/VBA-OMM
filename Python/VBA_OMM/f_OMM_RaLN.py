import numpy as np


def f_model(X, th, u, inF):
    # Function defining the ODE evolution of the van der Pol oscillator

    V = inF["V"]
    Gb = inF["Gb"]
    A = inF["A"]

    nf = 3
    nb = 5

    n_theta = nf + nb
    n = 2

    t = u[0, 0]
    I = u[1, 0]
    Rap = u[2, 0]

    Ra = f_Ra(t, th, A)[0]

    # Model Equations
    dx = np.zeros((n, 1))
    dx[0, 0] = -X[0]*X[1] - np.exp(th[0])*(X[0]-Gb) + (Ra+Rap)/V
    dx[1, 0] = -np.exp(th[1])*(X[1] - np.exp(th[2])*I)

    # Derivatives of the ODEs w.r.t
        # Model states
    dFdX = np.zeros((n, n))
    dFdX[0, 0] = -X[1] - np.exp(th[0])
    dFdX[0, 1] = -X[0]
    dFdX[1, 0] = 0
    dFdX[1, 1] = -np.exp(th[1])
        # Evolution parameters
    dFdTh = np.zeros((n, n_theta))
    dFdTh[0, 0] = -np.exp(th[0])*(X[0]-Gb)
    dFdTh[1, 0] = 0
    dFdTh[0, 1] = 0
    dFdTh[1, 1] = -np.exp(th[1])*(X[1] - np.exp(th[2])*I)
    dFdTh[0, 2] = 0
    dFdTh[1, 2] = np.exp(th[1]+th[2])*I

    dFdth = getDth(dFdTh, t, th, V, A)

    temp = dFdth.T
    return dx, dFdX, dFdth


def f_obs(X, phi, u, inG):
    # Observation of the glucose state

    n_phi = 0       # No of Observation Parameters
    nY = 1          # No of Observations
    n = 2           # Model order

    # Observation Equation
    gx = np.zeros((nY, 1))
    gx[0, 0] = X[0]

    # Derivatives of the Observation equation w.r.t
        # - Model states
    dGdX = np.zeros((nY, n))
    dGdX[0, 0] = 1
    dGdX[0, 1] = 0
        # - Observation Parameters
    dGdPhi = np.zeros((nY, n_phi))

    return gx, dGdX, dGdPhi


def f_Ra(t, th, A):

    nf = 2

    if t == 0:
        Ra = 0
        f1 = 0
        f2 = 0
    else:
        Rh = 1 / (1 + np.exp(-th[nf + 5]))
        f1 = (1 - Rh) * A / (t * np.sqrt(np.exp(th[nf + 2]) * np.pi)) * np.exp(-(np.log(t / np.exp(th[nf + 1])) - np.exp(th[nf + 2]) / 2) ** 2 / np.exp(th[nf + 2]))
        f2 = A * Rh / (t * np.sqrt(np.exp(th[nf + 4]) * np.pi)) * np.exp(-(np.log(t / np.exp(th[nf + 3])) - np.exp(th[nf + 4]) / 2) ** 2 / np.exp(th[nf + 4]))

        Ra = f1 + f2

    return Ra, f1, f2


def getDth(dFdTh, t, th, V, A):

    nf = 2

    if t == 0:
        dFdTh[0, nf+1] = 0
        dFdTh[0, nf+2] = 0
        dFdTh[0, nf+3] = 0
        dFdTh[0, nf+4] = 0
        dFdTh[0, nf+5] = 0
    else:
        Rh = 1 / (1 + np.exp(-th[nf + 5]))
        fu1 = (1 - Rh) * A / (t * np.sqrt(np.exp(th[nf + 2]) * np.pi)) * np.exp(
            -(np.log(t / np.exp(th[nf + 1])) - np.exp(th[nf + 2]) / 2) ** 2 / np.exp(th[nf + 2]))
        fu2 = A * Rh / (t * np.sqrt(np.exp(th[nf + 4]) * np.pi)) * np.exp(
            -(np.log(t / np.exp(th[nf + 3])) - np.exp(th[nf + 4]) / 2) ** 2 / np.exp(th[nf + 4]))

        dFdTh[0, nf + 1] = 2 * fu1 * np.exp(-th[nf+2]) * (-np.exp(th[nf+2]) / 2 + np.log(t / np.exp(th[nf+1]))) / V
        dFdTh[0, nf + 2] = fu1 * (-np.exp(th[nf+2])/2 + np.log(t/np.exp(th[nf+1]))+np.exp(-th[nf+2])*(-np.exp(th[nf+2])/2 + np.log(t/np.exp(th[nf+1])))**2)/V - fu1/2/V

        dFdTh[0, nf + 3] = 2 * fu2 * np.exp(-th[nf+4]) * (-np.exp(th[nf+4])/2 + np.log(t/np.exp(th[nf+3])))/V
        dFdTh[0, nf + 4] = fu2 * (-np.exp(th[nf+4])/2 + np.log(t/np.exp(th[nf+3]))+np.exp(-th[nf+4])*(-np.exp(th[nf+4])/2+np.log(t/np.exp(th[nf+3])))**2)/V - fu2/2/V

        dFdTh[0, nf + 5] = fu2 * np.exp(-th[nf+5])/(1+np.exp(-th[nf+5]))/V - A*Rh/(t*np.sqrt(np.exp(th[nf+2])*np.pi))*np.exp(-(np.log(t/np.exp(th[nf+1]))-np.exp(th[nf+2])/2)**2/np.exp(th[nf+2]))*np.exp(-th[nf+5])/(1+np.exp(-th[nf+5]))/V
    return dFdTh

