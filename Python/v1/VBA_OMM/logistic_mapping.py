import numpy as np


def f_logistic_mapping(m_i, s_i, mod):

    if mod == 1:
        m_o, s_o = Log2Norm(m_i, s_i)
    else:
        m_o, s_o = Norm2Log(m_i, s_i)

    return m_o, s_o


def Log2Norm(m_i, s_i):

    m_o = -np.log(1/m_i - 1)
    sig = np.arange(0.01, 10, 0.01)
    CV = np.zeros(np.size(sig))
    for i in range(0, np.size(sig)-1):
        _, CV[i] = Norm2Log(m_o, sig[i])

    tmp = np.abs(s_i-CV)
    mx = np.amax(tmp)
    idx = np.where(tmp == mx)
    s_o = sig[idx]

    return m_o, s_o


def Norm2Log(m_i, s_i):

    a = 3/np.pi**2

    m_o = f_logistic(m_i)
    mea = f_logistic(m_i/np.sqrt(1+a*np.sqrt(s_i)))
    vari = mea*(1-mea)*(1-1/np.sqrt(1+a*np.sqrt(s_i)))
    s_o = np.sqrt(vari)/mea*100

    return m_o, s_o


def f_logistic(x):

    return 1/(1+np.exp(-x))

