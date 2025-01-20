import numpy as np
import scipy.optimize as optimize
import scipy.linalg as la
import scipy.special as spec


def ToGamma(mu, sig):

    if mu < 5E-3:
        raise Exception('ERROR in Calculation of Gamma distribution: mu must be greater than 5e-3')
    if sig/mu < 5E-3:
        raise Exception('ERROR in Calculation of Gamma distribution: sig/mu must be greater than 5e-3')
    if sig/mu > 5:
        raise Exception('ERROR in Calculation of Gamma distribution: sig/mu must be smaller than 5')

    # Upper bound for a
    a0 = 1/8*(1+np.sqrt(49+mu**4/sig**4+50*mu**2/sig**2)+mu**2/sig**2)

    # Find a by minimizing the function D for 1<a<a0
    res = optimize.minimize_scalar(D_a, bounds=(1, a0), args=(mu, sig), method='bounded',  options={'xatol': 1e-04, 'maxiter': 500})
    a = res.x

    # Calculate b
    b = (mu/(np.exp(spec.loggamma(a-0.5)-spec.loggamma(a))))**2

    return a, b


def D_a(a, mu, sig):

    S = np.exp(2 * (spec.loggamma(a - 0.5) - spec.loggamma(a)))
    D = mu ** 2 / S - sig ** 2 / ((1 / (a - 1)) - S)
    D = np.log(D ** 2 + 1)

    return D


def FromGamma(a, b):

    if a <= 1:
        raise Exception('ERROR in Coversion of Gamma distribution: a must be greater than 1')
    else:

        mu = np.sqrt(b)*np.exp(spec.loggamma(a-0.5)-spec.loggamma(a))

        sig = np.sqrt(b/(a-1)-b*np.exp(2*(spec.loggamma(a-0.5)-spec.loggamma(a))))

        return mu, sig
