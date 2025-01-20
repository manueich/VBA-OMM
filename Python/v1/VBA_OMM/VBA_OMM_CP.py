import numpy as np
from . import VBA
from . import f_OMM_CP as f_CP
from . import logistic_mapping as LogMap
import matplotlib.pyplot as plt

def main(dat, priors, const, opt):

    # ----------------------------------------------
    # Check Inputs
    # Check dat
    try:
        t = np.reshape(dat["t"], (1, np.size(dat["t"])))
        G = np.reshape(dat["G"], (1, np.size(dat["G"])))
        CP = np.reshape(dat["CP"], (1, np.size(dat["CP"])))
    except:
        print("ERROR: Dat structure is flawed")
        return []
    if t[0, 0] != 0:
        print("ERROR: First datapoint must me at t=0")
        return []
    if t[0, 0] >= t[0, -1]:
        print("ERROR: The given time points are not ascending")
        return []

    # Check const
    try:
        age = const["age"]
        subject_type = const["subject_type"]
        dt = const["dt"]
        measCV = const["measCV"]
        CPb = const["CPb"]
        Gb = const["Gb"]
    except:
        print("ERROR: const structure is flawed")
        return []

    # Check Priors
    try:
        p_T = priors["T"]
        p_beta = priors["beta"]
        p_h = priors["h"]
        p_kd = priors["kd"]
    except:
        print('ERROR: Priors are flawed')
        return []

    # Check opt
    try:
        displayWin = opt["displayWin"]
        updateMeasCV = opt["updateMeasCV"]
    except:
        print("ERROR: opt structure is flawed")
        return []

    # --------------------------------------
    # Integration Time grid
    ti = np.arange(0, t[0, -1] + dt, dt)
    ti = ti.reshape(1, ti.size)

    # Input
    u = np.zeros((3, np.shape(ti)[1]))
    u[0, :] = ti
    Gi = np.interp(np.reshape(ti, np.shape(ti)[1]),
                   np.reshape(t, np.shape(t)[1]),
                   np.reshape(G, np.shape(G)[1]))
    u[1, :] = np.reshape(Gi, (1, np.shape(Gi)[0]))
    u[2, :] = np.diff(Gi,append=Gi[-1]) / dt

    data = {"y": CP[:, 0:np.size(CP)] - CPb,
            "t": t[:, 0:np.size(CP)],
            "u": u}

    k01,k12,k21 = get_CP_params(age,subject_type)

    inF = {"k01": k01,
           "k12": k12,
           "k21": k21}

    # Construct Priors
        # - Parameters
    muTheta = np.zeros((4,1))
    muTheta[0, 0] = np.log(p_T[0])
    muTheta[1, 0] = np.log(p_beta[0])
    muTheta[2, 0] = np.log(p_h[0])
    muTheta[3, 0] = np.log(p_kd[0])

    SigmaTheta = np.zeros((4,4))
    SigmaTheta[0, 0] = np.log((p_T[1] / 100) ** 2 + 1)
    SigmaTheta[1, 1] = np.log((p_beta[1] / 100) ** 2 + 1)
    SigmaTheta[2, 2] = np.log((p_h[1] / 100) ** 2 + 1)
    SigmaTheta[3, 3] = np.log((p_kd[1] / 100) ** 2 + 1)

        # - Initial Conditions
    muX0 = np.zeros((3, 1))
    muX0[0, 0] = 0
    muX0[1, 0] = 0
    muX0[2, 0] = 0

    SigmaX0 = np.zeros((3, 3))

        # - Measurement Noise
    CPi = np.interp(np.reshape(ti, np.shape(ti)[1]),
                   np.reshape(t, np.shape(t)[1]),
                   np.reshape(CP, np.shape(CP)[1]))

    iQy = [np.eye(1)*1/(CP[0, 1]/np.mean(CPi))]
    for i in range(0, np.size(t) - 1):
        iQy.append(np.eye(1)*1/(CP[0, i+1]/np.mean(CPi)))

    a, b = VBA.VBA_NoiseDistConv.ToGamma(np.mean(CPi) * measCV / 100, np.mean(CPi) * measCV / 100 * 0.5)

    pr = {"a": a,
          "b": b,
          "muTheta": muTheta,
          "SigmaTheta": SigmaTheta,
          "muX0": muX0,
          "SigmaX0": SigmaX0,
          "iQy": iQy}

    # Model Dimensions
    dim = {"n": 3,
           "n_theta": np.size(pr["muTheta"]),
           "n_phi": 0}

    options = {"f_model": f_CP.f_model,
               "f_obs": f_CP.f_obs,
               "inF": inF,
               "dim": dim,
               "verbose": False,
               "Display": False,
               "updateHP": updateMeasCV,
               "TolFun": 1E-4,
               "GnTolFun": 1e-4,
               "MaxIter": 100,
               "ODESolver": "RK"}

    # -----------------------------------------
    # ---- INVERSION --------

    if displayWin:
        print("Model Inversion ...")

    try:
        post, out_TB = VBA.VBA_Main.main(data, ti, pr, options)
    except:
        print("Inverison Failed")
        return []

    if displayWin:
       print("DONE")

    # -----------------------------------------
    # ---- WRAPUP RESULTS --------

    out = {"priors": priors,
           "options": opt,
           "const": const,
           "data": dat,
           "Toolbox": {"posterior": post, "out": out_TB}}

    # Posterior
    mu, sig = VBA.VBA_NoiseDistConv.FromGamma(post["a"],post["b"])    
    posterior = {"T": np.array([np.exp(post["muTheta"][0, 0]), np.sqrt(np.exp(post["SigmaTheta"][0, 0]) - 1) * 100]),
                 "beta": np.array([np.exp(post["muTheta"][1, 0]), np.sqrt(np.exp(post["SigmaTheta"][1, 1]) - 1) * 100]),
                 "h": np.array([np.exp(post["muTheta"][2, 0]), np.sqrt(np.exp(post["SigmaTheta"][2, 2]) - 1) * 100]),
                 "kd": np.array([np.exp(post["muTheta"][3, 0]), np.sqrt(np.exp(post["SigmaTheta"][3, 3]) - 1) * 100]),
                 "measCV": np.array([mu/np.mean(CPi)*100,sig/mu*100])}

    out.update({"posterior": posterior})

    # Model Output
    X, SigX, SR, SRs, SRd = get_ModelOut(post, out_TB, ti)

    Model_Output = {"t": ti,
                    "CP1": X[[0], :] + CPb,
                    "CP2": X[[1], :],
                    "SigCP1": SigX[[0], :],
                    "SigCP2": SigX[[1], :],
                    "SR": SR,
                    "SRs": SRs,
                    "SRd": SRd,                    
                    "PhiS": posterior["beta"][0],
                    "PhiD": posterior["kd"][0],
                    "PhiB": k01*CPb/Gb*1e3,
                    "Phi": posterior["beta"][0] + posterior["kd"][0]*(np.max(G) - Gb)/(np.trapz(Gi,ti)*dt)}
    out.update({"Model_Output": Model_Output})

    # Model performance
    CP_pred = out_TB["suffStat"]["gx"] + CPb
    CP_dat = data["y"] + CPb
    wres = (CP_dat - CP_pred) / (posterior["measCV"][0]*CP_dat/100)
    RMSE = (CP_dat - CP_pred) ** 2
    RMSE = np.sqrt(np.sum(RMSE) / np.size(RMSE))

    Performance = {"FreeEnergy": out_TB["F"],
                   "R2": out_TB["fit"]["R2"],
                   "AIC": out_TB["fit"]["AIC"],
                   "BIC": out_TB["fit"]["BIC"],
                   "LL": out_TB["fit"]["LL"],
                   "wres": wres,
                   "RMSE": RMSE}

    out.update({"Performance": Performance})

    if opt["displayWin"]:
        create_ResultsFigure(out)

    return out


def get_ModelOut(post, out_TB, t):

    X = out_TB["ModelOut"]["muX"]
    n = len(out_TB["ModelOut"]["SigmaX"])
    
    SigX = np.zeros((3, n))
    for i in range(0, n):
        SigX[0, i] = np.sqrt(out_TB["ModelOut"]["SigmaX"][i][0, 0])
        SigX[1, i] = np.sqrt(out_TB["ModelOut"]["SigmaX"][i][1, 1])
        SigX[2, i] = np.sqrt(out_TB["ModelOut"]["SigmaX"][i][2, 2])

    # SR
    u = out_TB["data"]["u"]
    SRs, SRd, SR = np.empty(np.shape(u)[1]), np.empty(np.shape(u)[1]), np.empty(np.shape(u)[1])
    for i in range(np.shape(u)[1]):
        SR[i], SRs[i], SRd[i] = f_CP.f_SR(X[:,i],post["muTheta"],u[2,i])

    return X, SigX, SR, SRs, SRd


def create_ResultsFigure(out):

    opt = out["options"]
    post = out["posterior"]
    pr = out["priors"]
    data = out["data"]
    mo = out["Model_Output"]    

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 6),tight_layout=True)
    
    # Summary
    ax[0, 0].set_title('SUMMARY')
    ax[0, 0].set_xticks([])
    ax[0, 0].set_yticks([])
    ax[0, 0].axis("off")

    strgs = ["Results of the C-Peptide OMM"]
    strgs.append("Model Fit criteria")
    strgs.append("   - R2: " + str(np.round(out["Performance"]["R2"], 3)))
    strgs.append("   - RMSE: " + str(np.round(out["Performance"]["RMSE"], 2)) + " mmol/L")
    strgs.append(r"Parameters [median $\pm$ CV [%]]")
    # strgs.append("   - T [min]: " + str(np.round(post["T"][0], 1)) + r" $\pm$ "
    #              + str(np.round(post["T"][1], 1)))
    # strgs.append("   - h [mmol/l]: " + str(np.round(post["h"][0], 1)) + r" $\pm$ "
    #              + str(np.round(post["h"][1], 1)))
    strgs.append(r"   - Static $\beta$-cell Senst.$\Phi_{s}$ [$10^{-9}$min$^{-1}$]:")
    strgs.append("       "+str(np.round(post["beta"][0], 1)) + r" $\pm$ " + str(np.round(post["beta"][1], 1)))
    strgs.append(r"   - Dynamic $\beta$-cell Senst.$\Phi_{d}$ [$10^{-9}$]:")
    strgs.append("       "+str(np.round(post["kd"][0], 1)) + r" $\pm$ " + str(np.round(post["kd"][1], 1)))
    strgs.append(r"   - Basal $\beta$-cell Senst.$\Phi_{b}$ [$10^{-9}$min$^{-1}$]:")
    strgs.append("       "+str(np.round(mo["PhiB"],1)))
    strgs.append(r"   - Total $\beta$-cell Senst.$\Phi$ [$10^{-9}$min$^{-1}$]:")
    strgs.append("       "+str(np.round(float(mo["Phi"]), 1)))


    for i in range(0, len(strgs)):
        ax[0, 0].text(0.05, 0.9-i/13, strgs[i], color='k')

    # Data vs Model
    ax[0, 1].set_title('Data vs Model Prediction')
    ax[0, 1].set_ylabel("C-Peptide [nmol/L]")
    ax[0, 1].set_xlabel("Time [min]")
    ax[0, 1].axis()

    # Data
    h1 = ax[0, 1].plot(To1D(data["t"]), To1D(data["CP"]), marker='o', markersize=4, ls='--', linewidth=0.5, color="black",label="Data")
    # Model Pred
    h3 = ax[0, 1].plot(To1D(mo["t"]), To1D(mo["CP1"]), color="tab:blue")
    h4 = ax[0, 1].fill_between(To1D(mo["t"]), To1D(mo["CP1"] - mo["SigCP1"]), To1D(mo["CP1"] + mo["SigCP1"]), alpha=0.2, color="tab:blue")
    ax[0,1].legend(handles=[h1[0],(h3[0],h4)],labels=["Data",r"Model Pred. $\pm1\sigma$"])

    # Secreation Rate
    ax[1, 0].set_title('Secretion rate')
    ax[1, 0].set_ylabel("SR [nmol/l/min]")
    ax[1, 0].set_xlabel("Time [min]")
    ax[1, 0].axis()

    ax[1, 0].plot(To1D(mo["t"]), To1D(mo["SRd"]), color="tab:red",label="SRd")
    ax[1, 0].plot(To1D(mo["t"]), To1D(mo["SRs"]), color="tab:green",label="SRs")
    ax[1, 0].plot(To1D(mo["t"]), To1D(mo["SR"]), color="tab:blue",label="SR")
    ax[1, 0].legend() 

    # Weighted Residuals
    ax[1, 1].set_title('Weighted Residuals')
    ax[1, 1].set_ylabel("WRES")
    ax[1, 1].set_xlabel("Time [min]")
    ax[1, 1].axis()

    ax[1, 1].plot(To1D(data["t"]), 0 * np.ones(np.size(data["t"])), ls='-', linewidth=1,
                  color='k')
    ax[1, 1].plot(To1D(data["t"]), -1 * np.ones(np.size(data["t"])), ls='--', linewidth=1,
                  color='k')
    ax[1, 1].plot(To1D(data["t"]), 1 * np.ones(np.size(data["t"])), ls='--', linewidth=1,
                  color='k')
    ax[1, 1].plot(To1D(data["t"]), To1D(out["Performance"]["wres"]), marker='o', markersize=3, ls='--', linewidth=1,
                  color="tab:blue")

    plt.show()


def get_CP_params(age,subject_type):
    """
    The method to calculate C-Petide parameters is described in
    van Cauter et al., DIABETES, VOL 41, MARCH 1992
    """

    if subject_type == "normal":
        HLs = 4.95
        f = 0.76
    elif subject_type == "obese":
        HLs = 4.55
        f = 0.78
    elif subject_type == "NIDDM":
        HLs = 4.52
        f = 0.78
    HLl = 0.14 * age + 29.2

    a = np.log(2)/HLs
    b = np.log(2)/HLl

    k12 = (b+(1/f-1)*a)*f
    k01 = a*b/k12
    k21 = a+b-k12-k01

    return k01,k12,k21


def To1D(x):

    return np.reshape(x, (np.size(x)))


def getLNbounds (x, sc):

    mu = float(x[0]) * sc
    sig = float(np.exp(np.sqrt(np.log((x[1]/100)**2 + 1))))

    up = mu*sig
    lo = mu/sig
    dat = np.array([[mu-lo], [up-mu]])

    return dat

def getLogisticbounds (x):

    m_o, s_o = LogMap.f_logistic_mapping(x[0], x[1], 1)
    lo = LogMap.f_logistic(m_o - np.sqrt(s_o))
    up = LogMap.f_logistic(m_o + np.sqrt(s_o))
    mu = x[0]

    dat = np.array([[float(mu-lo)], [float(up-mu)]])

    return dat