import numpy as np


def f_model(X, th, u, inF):
    # Function defining the ODE evolution of the van der Pol oscillator

    V = inF["V"]
    Gb = inF["Gb"]
    A = inF["A"]
    al = inF["alpha"]
    tb = inF["tb"]

    nf = 3
    nb = np.shape(tb)[0]-2

    n_theta = nf + nb
    n = 2

    t = u[0, 0]
    I = u[1, 0]
    Rap = u[2, 0]

    Ra = f_Ra(t, tb, th, A, al)

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

    dFdth = getDth(dFdTh, t, th, tb, al, V)

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


def f_Ra(t, tb, th, A, al):

    k = np.exp(th[3:np.shape(th)[0], 0])
    s = (-2 * A * al + al * (-k[1] * tb[1] + k[0] * tb[2] - k[2] * tb[2] + k[1] * tb[3] - k[3] * tb[3] + k[2] * tb[4]
                             - k[4] * tb[4] + k[3] * tb[5] + k[4] * tb[6]) + k[5] * (2 - al * tb[6] + al * tb[7])) / (al * (tb[5] - tb[7]))

    k = np.concatenate((np.array([0]), k[0:5], np.array([s]), k[5:np.shape(th)[0]]))

    for i in range(1, np.shape(tb)[0]):
        if tb[i - 1] <= t <= tb[i]:
            Ra = k[i-1] + (k[i] - k[i-1]) / (tb[i] - tb[i-1]) * (t - tb[i-1])

    if t > tb[np.shape(tb)[0]-1]:
        Ra = k[np.shape(k)[0]-1] * np.exp(-(t - tb[np.shape(tb)[0]-1]) * al)

    return Ra


def getDth(dFdTh, t, th, tb, al, V):

    nf = 3

    k = np.exp(th[3:np.shape(th)[0], 0])

    if tb[0] <= t <= tb[1]:
        dFdTh[0, nf] = k[0] * t / tb[1] / V

    if tb[1] <= t <= tb[2]:
        dFdTh[0, nf] = k[0] * (tb[2] - t) / (tb[2] - tb[1]) / V
        dFdTh[0, nf+1] = k[1] * (t - tb[1]) / (tb[2] - tb[1]) / V

    if tb[2] <= t <= tb[3]:
        dFdTh[0, nf+1] = k[1] * (tb[3] - t) / (tb[3] - tb[2]) / V
        dFdTh[0, nf+2] = k[2] * (t - tb[2]) / (tb[3] - tb[2]) / V

    if tb[3] <= t <= tb[4]:
        dFdTh[0, nf+2] = k[2] * (tb[4] - t) / (tb[4] - tb[3]) / V
        dFdTh[0, nf+3] = k[3] * (t - tb[3]) / (tb[4] - tb[3]) / V

    if tb[4] <= t <= tb[5]:
        dFdTh[0, nf+3] = k[3] * (tb[5] - t) / (tb[5] - tb[4]) / V
        dFdTh[0, nf+4] = k[4] * (t - tb[4]) / (tb[5] - tb[4]) / V

    if tb[5] <= t <= tb[6]:
        dFdTh[0, nf] = k[0] * tb[2] * (t - tb[5]) / (tb[6] - tb[5]) / (tb[5] - tb[7]) / V
        dFdTh[0, nf+1] = k[1] * (tb[1] - tb[3]) * (t - tb[5]) / (tb[5] - tb[6]) / (tb[5] - tb[7]) / V
        dFdTh[0, nf+2] = k[2] * (tb[2] - tb[4]) * (t - tb[5]) / (tb[5] - tb[6]) / (tb[5] - tb[7]) / V
        dFdTh[0, nf+3] = -k[3] * (tb[5] - tb[3]) * (t - tb[5]) / (tb[5] - tb[6]) / (tb[5] - tb[7]) / V
        dFdTh[0, nf+4] = k[4] * (-tb[4] * tb[5] + t * (tb[4] + tb[5] - tb[6] - tb[7]) + tb[6] * tb[7]) / (tb[5] - tb[6]) / (tb[5] - tb[7]) / V
        dFdTh[0, nf+5] = k[5] * (t - tb[5]) * (2 + al * (tb[7] - tb[6])) / (tb[6] - tb[5]) / (tb[5] - tb[7]) / al / V

    if tb[6] <= t <= tb[7]:
        dFdTh[0, nf] = k[0] * tb[2] * (tb[7] - t) / (tb[7] - tb[6]) / (tb[5] - tb[7]) / V
        dFdTh[0, nf+1] = k[1] * (tb[1] - tb[3]) * (t - tb[7]) / (tb[5] - tb[7]) / (tb[7] - tb[6]) / V
        dFdTh[0, nf+2] = k[2] * (tb[2] - tb[4]) * (t - tb[7]) / (tb[5] - tb[7]) / (tb[7] - tb[6]) / V
        dFdTh[0, nf+3] = k[3] * (tb[3] - tb[5]) * (t - tb[7]) / (tb[5] - tb[7]) / (tb[7] - tb[6]) / V
        dFdTh[0, nf+4] = k[4] * (tb[4] - tb[6]) * (t - tb[7]) / (tb[6] - tb[7]) / (tb[7] - tb[5]) / V
        dFdTh[0, nf+5] = k[5] * (t * (-2 + al * (tb[5] + tb[6] - 2 * tb[7])) + 2 * tb[7] + al * (-tb[5] * tb[6] + tb[7]**2)) / (tb[5] - tb[7]) / (tb[7] - tb[6]) / al / V

    return dFdTh

