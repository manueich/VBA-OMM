"""
VBA_OMM:
A Python package for the inversion of the oral minimal model (OMM) from non-fasting conditions
using variational Bayesian analysis

-----------------------------------------------
Glucose OMM
This package inverts the OMM of glucose dynamics in the following form
      dG/dt = - G*X - p1*(G-Gb) + (Ra+Rap)/V      G(0) = G0
      dX/dt = - p2* (X - SI*(I-Ib))               X(0) = X0

The main function of the package is called as follows:

         out = VBA_OMM.mainG(dat, priors, const, opt)

INPUT
  - dat: Dictionary containing the data
      - t: ndarray of Time of sampling points in min. The first datapoint must be at t=0
      - G: ndarray of Glucose data at time points t in mmol/L
      - I: ndarray of Insulin data at time points t in arbitrary insulin unit (IU),
      e.g. mU/L, or pmol/L
  - opt: Dictionary specifying inversion options
      - GA_fun: Either 'RaPL5', 'RaPL7', or 'RaLN', using either the piecewise
      linear or log-normally based functions representing Ra
      - tb: ndarray of Time of breakpoints of RaPL in min. First breakpoint must be 0.
      - alpha: Exponential decay rate of RaPL after last breakpoint in
      1/min
      - displayWin: True or False specifying whether inversion results are
      displayed
  - priors: Dictionary specifying the priors. A coefficient of variation
  (CV) of zero means the parameter is fixed and not updated during
  inversion
      - p1: 1x2 ndarray specifying median and CV of log-normal distribution
      in 1/min and %
      - p2: 1x2 ndarray specifying median and CV of log-normal distribution
      in 1/min and %
      - SI: 1x2 ndarray specifying median and CV of log-normal distribution
      in 1/min per IU and %
      - k: Mx2 ndarray specifying median and CV of Ra function parameters.
      CV in %.
      For RaPL5: M=4. Heights of Ra at breakpoints tb in
      mmol/kg/min, starting  with tb(2). The height at the penultimate
      breakpoint is calculated from the total AUC of RA and therefore
      not specified.     
      For RaPL7: M=6. Heights of Ra at breakpoints tb in
      mmol/kg/min, starting  with tb(2). The height at the penultimate
      breakpoint is calculated from the total AUC of RA and therefore
      not specified.       
      For RaLN: M=5. Representing in order T1 in min, W1 no unit,
      T2 in min, W2 no unit and RH no unit. RH is restricted to (0,1)
      and is not log-normally distributed.
  - const: Dictionary specifying model constants
      - dt: ODE integration step size in min
      - A: AUC of Ra in mmol/kg, typically calculated as from the CHO content of the meal times 0.9
      - V: Glucose distribution volume in L/kg
      - X0: Initial condition of X, i.e. X(0) in 1/min
      - G0: Initial condition of G, i.e. G(0) in mmol/L
      - Gb: Basal level of G in mmol/L
      - Ib: Basal level of I in IU
      - measCV: measurement uncertainty CV of glucose assay in %
      - Rap: Persisting absorption from a previous meal. An empty input
      means no persisting absoption. If present, Rap has to be a ndarray
      coinciding with the integration time points on the
      grid t(1):dt:t(end)

OUTPUT out: Dictionary containing inversion results and input
  - priors, options, data, const: see INPUT
  - posterior: Dictionary containing posterior parameter distributions
  analogous to prior dictionary
  - Model_Output: Dictionary containing model output
      - t: ODE integration time grid in min
      - G, SigG: Model inferred glucose and uncertainty (SD) on t
      in mmol/L
      - X, SigX: Model inferred state X and uncertainty (SD) on t
      in 1/min
      - Ra, SigRa: Model inferred glucose appearance and uncertainty (SD)
      on t in in mmol/kg/min
      - Rap: Persisting appearance for possible consecutive meal from
      t(end) to 2*t(end)
  - Performance: Dictionary containing model performance metrics
      - FreeEnergy: Variational free energy, i.e. lower bound on log
      model evidence
      - R2: Coefficient of determination between data and model output
      - AIC: Akaike Information Criterion
      - BIC: Bayesian Information Criterion
      - LL: Log Likelihood
      - RMSE: Root mean squared error between data and model output
      - wres: Weighted residuals between data and model output
  - VB_Toolbox: Raw output of the VB method. See https://github.com/manueich/VBA-python for details and documentation

------------------------------------------------
C-Peptide OMM
This package inverts the OMM of C-peptide dynamics in the following form
      dCP1/dt = -(k01+k21)*CP1 + k12*CP2 + SR     CP1(0) = 0   
      dCP2/dt = k21*CP1 -k12*CP2                  CP2(0) = 0
      SR = SRs + SRd
      dSRs/dt = -1/T*[SRs - beta*(G - h)]
      SRd = kd*dG/dt if G > 0
      SRd = 0        if G <= 0

A detailed description of the model is provided in 
Sunehag et al., Obesity (2008) 17, 233–239. doi:10.1038/oby.2008.496

The main function of the package is called as follows:

         out = VBA_OMM.mainCP(dat, priors, const, opt)

INPUT
  - dat: Dictionary containing the data
      - t: ndarray of Time of sampling points in min. The first datapoint must be at t=0
      - G: ndarray of Glucose data at time points t in mmol/L
      - CP: ndarray of C-Peptide data at time points t in nmol/L
  - opt: Dictionary specifying inversion options
      - displayWin: True or False specifying whether inversion results are
      displayed
      - updateMeasCV: True or False whether MeasCV is updated during inversion. If so,
      the results are stored in the posterior dictionary
  - priors: Dictionary specifying the priors. A coefficient of variation
  (CV) of zero means the parameter is fixed and not updated during
  inversion
      - T: 1x2 ndarray specifying median and CV of log-normal distribution
      in min and %
      - beta: 1x2 ndarray specifying median and CV of log-normal distribution
      in 1E-9 1/min and %
      - h: 1x2 ndarray specifying median and CV of log-normal distribution
      in mmol/L and %
      - kd: 1x2 ndarray specifying median and CV of log-normal distribution
      in 1E-9 and %
  - const: Dictionary specifying model constants
      - dt: ODE integration step size in min
      - Gb: Basal level of G in mmol/L
      - CPb: Basal level of CP in nmol/L
      - measCV: measurement uncertainty CV of C-Peptide assay in %
      - age: Subject age to calculate dynamic C-petide parameters (k01,k12,k21)
      - subject_type: Glycemic status of the subjects "normal", "weight" or "NIDDM"
      The method to calculate C-Petide parameters is described in
      van Cauter et al., DIABETES, VOL 41, MARCH 1992

OUTPUT out: Dictionary containing inversion results and input
  - priors, options, data, const: see INPUT
  - posterior: Dictionary containing posterior parameter distributions
  analogous to prior dictionary, except for measCV which is provided as 
  CV [%] and CV [%] (uncertainty in measCV)
  - Model_Output: Dictionary containing model output
      - t: ODE integration time grid in min
      - CP1, SigCP1: Model inferred c-peptide1 and uncertainty (SD) on t
      in nmol/L
      - CP2, SigCP2: Model inferred c-peptide2 and uncertainty (SD) on t
      in nmol/L
      - SR: Total secretion rate in nmol/kg/min
      - SRs: Static secretion rate in nmol/kg/min
      - SRd: Dynamic secretion rate in nmol/kg/min
      - PhiS: Static beta cell responsivity (beta) in 1E-9 1/min
      - PhiD: Dynamic beta cell responsivity (kd) in 1E-9
      - PhiB: Basal beta cell responsivity in 1E-9 1/min
      - Phi: Total beta cell responsivity in 1E-9 1/min
  - Performance: Dictionary containing model performance metrics
      - FreeEnergy: Variational free energy, i.e. lower bound on log
      model evidence
      - R2: Coefficient of determination between data and model output
      - AIC: Akaike Information Criterion
      - BIC: Bayesian Information Criterion
      - LL: Log Likelihood
      - RMSE: Root mean squared error between data and model output
      - wres: Weighted residuals between data and model output
  - VB_Toolbox: Raw output of the VB method. See https://github.com/manueich/VBA-python for details and documentation

-------------------------

Data1.csv
Data taken from figures in Breda et al., DIABETES, VOL. 50, JANUARY 2001
Data2.csv
Data taken from figures in Sunehag et al., Obesity (2008) 17, 233–239. doi:10.1038/oby.2008.496

M. Eichenlaub 16/10/2022
"""

from . import VBA_OMM_CP
from . import VBA_OMM_G,VBA_OMM_CP


def mainG(dat, priors, const, opt):

    return VBA_OMM_G.main(dat, priors, const, opt)

def mainCP(dat, priors, const, opt):

    return VBA_OMM_CP.main(dat, priors, const, opt)