# VBA_OMM:
# A Python package for the inversion of the oral minimal model (OMM) of glucose dynamics from non-fasting conditions
# using variational Bayesian analysis
#
# This package inverts the OMM in the following form
#       dG/dt = - G*X - p1*(G-Gb) + (Ra+Rap)/V      G(0) = G0
#       dX/dt = - p2* (X - SI*(I-Ib))               X(0) = X0
#
# The main function of the package is called as follows:
#
#          out = VBA_OMM.main(dat, priors, const, opt)
#
# INPUT
#   - dat: Dictionary containing the data
#       - t: ndarray of Time of sampling points in min. The first datapoint must be at t=0
#       - G: ndarray of Glucose data at time points t in mmol/L
#       - I: ndarray of Insulin data at time points t in arbitrary insulin unit (IU),
#       e.g. mU/L
#   - opt: Dictionary specifying inversion options
#       - GA_fun: Either 'RaPL' or 'RaLN', using either the piecewise
#       linear or log-normally based functions representing Ra
#       - tb: ndarray of Time of breakpoints of RaPL in min. First breakpoint must be
#       0 and last breakpoint must coincide with last datapoint.
#       - alpha: Exponential decay rate of RaPL after last breakpoint in
#       1/min
#       - displayWin: True or False specifying whether inversion results are
#       displayed
#   - priors: Dictionary specifying the priors. A coefficient of variation
#   (CV) of zero means the parameter is fixed and not updated during
#   inversion
#       - p1: 1x2 ndarray specifying median and CV of log-normal distribution
#       in 1/min and %
#       - p2: 1x2 ndarray specifying median and CV of log-normal distribution
#       in 1/min and %
#       - SI: 1x2 ndarray specifying median and CV of log-normal distribution
#       in 1/min per IU and %
#       - k: Mx2 ndarray specifying median and CV of Ra function parameters.
#       CV in %.
#       For RaPL: M=6. Heights of Ra at breakpoints tb in
#       mmol/kg/min, starting  with tb(2). The height at the penultimate
#       breakpoint is calculated from the total AUC of RA and therefore
#       not specified. NB: the package currently only supports RaPL with
#       exactly 8 breakpoints. The times of these breakpoints can however
#       be chosen freely. Please contact the developers if a different
#       number of breakpoints is required.
#       For RaLN: M=5. Representing in order T1 in min, W1 no unit,
#       T2 in min, W2 no unit and RH no unit. RH is restricted to (0,1)
#       and is not log-normally distributed.
#   - const: Dictionary specifying model constants
#       - dt: ODE integration step size in min
#       - A: AUC of Ra in mmol/kg
#       - V: Glucose distribution volume in L/kg
#       - X0: Initial condition of X, i.e. X(0) in 1/min
#       - G0: Initial condition of G, i.e. G(0) in mmol/L
#       - Gb: Basal level of G in mmol/L
#       - Ib: Basal level of I in IU
#       - measCV: measurement uncertainty CV of glucose assay in %
#       - Rap: Persisting absorption from a previous meal. An empty input
#       means no persisting absoption. If present, Rap has to be a ndarray
#       coinciding with the integration time points on the
#       grid t(1):dt:t(end)
#
# OUTPUT out: Dictionary containing inversion results and input
#   - priors, options, data, const: see INPUT
#   - posterior: Dictionary containing posterior parameter distributions
#   analogous to prior dictionary
#   - Model_Output: Dictionary containing model output
#       - t: ODE integration time grid in min
#       - G, SigG: Model inferred glucose and uncertainty (SD) on t
#       in mmol/L
#       - X, SigX: Model inferred state X and uncertainty (SD) on t
#       in 1/min
#       - Ra, SigRa: Model inferred glucose appearance and uncertainty (SD)
#       on t in in mmol/kg/min
#       - Rap: Persisting appearance for possible consecutive meal from
#       t(end) to 2*t(end)
#   - Performance: Dictionary containing model performance metrics
#       - FreeEnergy: Variational free energy, i.e. lower bound on log
#       model evidence
#       - R2: Coefficient of determination between data and model output
#       - AIC: Akaike Information Criterion
#       - BIC: Bayesian Information Criterion
#       - LL: Log Likelihood
#       - RMSE: Root mean squared error between data and model output
#       - wres: Weighted residuals between data and model output
#   - VB_Toolbox: Raw output of the VB method. See https://github.com/manueich/VBA-python for details and documentation
#
#   M. Eichenlaub 08/01/2020

from . import VBA_OMM_G


def main(dat, priors, const, opt):

    return VBA_OMM_G.main(dat, priors, const, opt)
