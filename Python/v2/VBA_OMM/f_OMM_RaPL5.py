import numpy as np

def exp(var):
    return np.exp(var)

def f_model(X, th, u, inF):
    # Function defining the OMM using RaPL

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

    s = (-2*A*al - al*tb[0]*exp(th[3]) - al*tb[1]*exp(th[4]) + al*tb[2]*exp(th[3]) - al*tb[2]*exp(th[5]) + al*tb[3]*exp(th[4]) + al*tb[4]*exp(th[5]) - al*tb[4]*exp(th[6]) + al*tb[5]*exp(th[6]) + 2*exp(th[6]))/(al*(tb[3] - tb[5]))


    if tb[0] <= t <= tb[1]:
        Ra = exp(th[3]) * (t - tb[0]) / (tb[1] - tb[0])
    elif tb[1] <= t <= tb[2]:
        Ra = exp(th[3]) + (exp(th[4]) - exp(th[3])) * (t - tb[1]) / (tb[2] - tb[1])
    elif tb[2] <= t <= tb[3]:
        Ra = exp(th[4]) + (exp(th[5]) - exp(th[4])) * (t - tb[2]) / (tb[3] - tb[2])
    elif tb[3] <= t <= tb[4]:
        Ra = exp(th[5]) + (s - exp(th[5])) * (t - tb[3]) / (tb[4] - tb[3])
    elif tb[4] <= t <= tb[5]:
        Ra = s + (exp(th[6]) - s) * (t - tb[4]) / (tb[5] - tb[4])
    else:
        Ra = exp(th[6]) * np.exp(-al*(t-tb[5]))

    return Ra

def getDth(dFdTh, t, th, tb, al, V):

    nf = 3   
    if tb[0] <= t <= tb[1]:
        dFdTh[0, nf] = (t - tb[0])*exp(th[3])/(V*(-tb[0] + tb[1]))
    elif tb[1] <= t <= tb[2]:
        dFdTh[0, nf] = (-(t - tb[1])*exp(th[3])/(-tb[1] + tb[2]) + exp(th[3]))/V
        dFdTh[0, nf+1] = (t - tb[1])*exp(th[4])/(V*(-tb[1] + tb[2]))
    elif tb[2] <= t <= tb[3]:
        dFdTh[0, nf+1] = (-(t - tb[2])*exp(th[4])/(-tb[2] + tb[3]) + exp(th[4]))/V
        dFdTh[0, nf+2] = (t - tb[2])*exp(th[5])/(V*(-tb[2] + tb[3]))
    elif tb[3] <= t <= tb[4]:
        dFdTh[0, nf] = (t - tb[3])*(-al*tb[0]*exp(th[3]) + al*tb[2]*exp(th[3]))/(V*al*(-tb[3] + tb[4])*(tb[3] - tb[5]))
        dFdTh[0, nf+1] = (t - tb[3])*(-al*tb[1]*exp(th[4]) + al*tb[3]*exp(th[4]))/(V*al*(-tb[3] + tb[4])*(tb[3] - tb[5]))
        dFdTh[0, nf+2] = ((t - tb[3])*(-exp(th[5]) + (-al*tb[2]*exp(th[5]) + al*tb[4]*exp(th[5]))/(al*(tb[3] - tb[5])))/(-tb[3] + tb[4]) + exp(th[5]))/V
        dFdTh[0, nf+3] = (t - tb[3])*(-al*tb[4]*exp(th[6]) + al*tb[5]*exp(th[6]) + 2*exp(th[6]))/(V*al*(-tb[3] + tb[4])*(tb[3] - tb[5]))
    elif tb[4] <= t <= tb[5]:
        dFdTh[0, nf] = (-(t - tb[4])*(-al*tb[0]*exp(th[3]) + al*tb[2]*exp(th[3]))/(al*(tb[3] - tb[5])*(-tb[4] + tb[5])) + (-al*tb[0]*exp(th[3]) + al*tb[2]*exp(th[3]))/(al*(tb[3] - tb[5])))/V
        dFdTh[0, nf+1] = (-(t - tb[4])*(-al*tb[1]*exp(th[4]) + al*tb[3]*exp(th[4]))/(al*(tb[3] - tb[5])*(-tb[4] + tb[5])) + (-al*tb[1]*exp(th[4]) + al*tb[3]*exp(th[4]))/(al*(tb[3] - tb[5])))/V
        dFdTh[0, nf+2] = (-(t - tb[4])*(-al*tb[2]*exp(th[5]) + al*tb[4]*exp(th[5]))/(al*(tb[3] - tb[5])*(-tb[4] + tb[5])) + (-al*tb[2]*exp(th[5]) + al*tb[4]*exp(th[5]))/(al*(tb[3] - tb[5])))/V
        dFdTh[0, nf+3] = ((t - tb[4])*(exp(th[6]) - (-al*tb[4]*exp(th[6]) + al*tb[5]*exp(th[6]) + 2*exp(th[6]))/(al*(tb[3] - tb[5])))/(-tb[4] + tb[5]) + (-al*tb[4]*exp(th[6]) + al*tb[5]*exp(th[6]) + 2*exp(th[6]))/(al*(tb[3] - tb[5])))/V
    else:
        dFdTh[0, nf+3] = exp(th[6])*exp(-al*(t - tb[5]))/V
    
    return dFdTh