import numpy as np
from . import VBA
from . import f_OMM_RaPL as f_RaPL
from . import f_OMM_RaLN as f_RaLN
from . import logistic_mapping as LogMap
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

def main(dat, priors, const, opt):

    # ----------------------------------------------
    # Check Inputs
    # Check dat
    try:
        t = np.reshape(dat["t"], (1, np.size(dat["t"])))
        G = np.reshape(dat["G"], (1, np.size(dat["G"])))
        I = np.reshape(dat["I"], (1, np.size(dat["I"])))
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
        A = const["A"]
        V = const["V"]
        dt = const["dt"]
        if not const["Rap"]:
            Rap = np.zeros((1, int(t[0, -1]/dt+1)))
        else:
            Rap = np.reshape(const["Rap"], (1, const["Rap"]))
            if np.size(Rap) != t[-1]/dt+1:
                print("ERROR: Rap is the wrong length")
                return []

        X0 = const["X0"]
        measCV = const["measCV"]
        G0 = const["G0"]
        Gb = const["Gb"]
        Ib = const["Ib"]
    except:
        print("ERROR: const structure is flawed")
        return []

    # Check Priors
    try:
        p_p1 = priors["p1"]
        p_p2 = priors["p2"]
        p_SI = priors["SI"]
        p_k = priors["k"]
    except:
        print('ERROR: Priors are flawed')
        return []

    # Check opt
    try:
        displayWin = opt["displayWin"]
        GA_fun = opt["GA_fun"]

        if GA_fun == 'RaPL':
            fname = f_RaPL.f_model
            gname = f_RaPL.f_obs
            tb = opt["tb"]
            alpha = opt["alpha"]
            if tb[0] != t[0, 0]:
                print("ERROR: First breakpoint of RaPL must coincide with first datapoint, i.e t=0")
                return []
            if tb[-1] != t[0, -1]:
                print("ERROR: Last breakpoint of RaPL must coincide with last datapoint")
                return []
            if np.shape(p_k)[0]+2 != np.size(tb):
                print("ERROR: Number of Breakpoints in RaPL does not match number of prior parameters")
                return []
        elif GA_fun == 'RaLN':
            fname = f_RaLN.f_model
            gname = f_RaLN.f_obs
            if np.shape(p_k)[0] != 5:
                print("ERROR: Wrong number of prior parameters for RaLN")
                return []
        else:
            print("ERROR: GA function type flawed")
            return []
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
    Ii = np.interp(np.reshape(ti, np.shape(ti)[1]),
                   np.reshape(t, np.shape(t)[1]),
                   np.reshape(I, np.shape(I)[1]))
    u[1, :] = np.reshape(Ii, (1, np.shape(Ii)[0])) - Ib
    u[2, :] = Rap

    data = {"y": G[:, 0:np.size(G)],
            "t": t[:, 0:np.size(G)],
            "u": u}

    inF = {"V": V,
           "Gb": Gb,
           "A": A}

    # Construct Priors
        # - Parameters
    muTheta = np.zeros((3 + np.shape(p_k)[0], 1))
    muTheta[0, 0] = np.log(p_p1[0])
    muTheta[1, 0] = np.log(p_p2[0])
    muTheta[2, 0] = np.log(p_SI[0])

    SigmaTheta = np.zeros((3 + np.shape(p_k)[0], 3 + np.shape(p_k)[0]))
    SigmaTheta[0, 0] = np.log((p_p1[1] / 100) ** 2 + 1)
    SigmaTheta[1, 1] = np.log((p_p2[1] / 100) ** 2 + 1)
    SigmaTheta[2, 2] = np.log((p_SI[1] / 100) ** 2 + 1)

    if GA_fun == 'RaPL':
        for i in range(0, np.shape(p_k)[0]):
            muTheta[3+i, 0] = np.log(p_k[i, 0])
            SigmaTheta[3+i, 3+i] = np.log((p_k[i, 1] / 100) ** 2 + 1)
            inF.update({"alpha": alpha,
                        "tb": tb})
    if GA_fun == 'RaLN':
        for i in range(0, np.shape(p_k)[0] - 1):
            muTheta[3+i, 0] = np.log(p_k[i, 0])
            SigmaTheta[3+i, 3+i] = np.log((p_k[i, 1] / 100) ** 2 + 1)
        muTheta[7, 0], SigmaTheta[7, 7] = LogMap.f_logistic_mapping(p_k[4, 0], p_k[4, 1], 1)

        # - Initial Conditions
    muX0 = np.zeros((2, 1))
    muX0[0, 0] = G0
    muX0[1, 0] = X0

    SigmaX0 = np.zeros((2, 2))

        # - Measurement Noise
    Gi = np.interp(np.reshape(ti, np.shape(ti)[1]),
                   np.reshape(t, np.shape(t)[1]),
                   np.reshape(G, np.shape(G)[1]))

    iQy = [np.eye(1)*1/(G[0, 1]/np.mean(Gi))]
    for i in range(0, np.size(t) - 1):
        iQy.append(np.eye(1)*1/(G[0, i+1]/np.mean(Gi)))

    a, b = VBA.VBA_NoiseDistConv.ToGamma(np.mean(Gi) * measCV / 100, np.mean(Gi) * measCV / 100 * 0.1)

    pr = {"a": a,
          "b": b,
          "muTheta": muTheta,
          "SigmaTheta": SigmaTheta,
          "muX0": muX0,
          "SigmaX0": SigmaX0,
          "iQy": iQy}

    # Model Dimensions
    dim = {"n": 2,
           "n_theta": np.size(pr["muTheta"]),
           "n_phi": 0}

    options = {"f_model": fname,
               "f_obs": gname,
               "inF": inF,
               "dim": dim,
               "verbose": False,
               "Display": False,
               "updateHP": False,
               "TolFun": 1E-4,
               "GnTolFun": 1e-4,
               "MaxIter": 100,
               "ODESolver": "Euler"}

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

    posterior = {"p1": np.array([np.exp(post["muTheta"][0, 0]), np.sqrt(np.exp(post["SigmaTheta"][0, 0]) - 1) * 100]),
                 "p2": np.array([np.exp(post["muTheta"][1, 0]), np.sqrt(np.exp(post["SigmaTheta"][1, 1]) - 1) * 100]),
                 "SI": np.array([np.exp(post["muTheta"][2, 0]), np.sqrt(np.exp(post["SigmaTheta"][2, 2]) - 1) * 100])}

    if GA_fun == 'RaPL':
        k = np.zeros((np.size(post["muTheta"]) - 3, 2))
        for i in range(0, np.size(post["muTheta"]) - 3):
            k[i, 0] = np.exp(post["muTheta"][i + 3])
            k[i, 1] = np.sqrt(np.exp(post["SigmaTheta"][i + 3, i + 3]) - 1) * 100
        # Correlation missing
    if GA_fun == 'RaLN':
        k = np.zeros((5, 2))
        for i in range(0, 3):
            k[i, 0] = np.exp(post["muTheta"][i + 3])
            k[i, 1] = np.sqrt(np.exp(post["SigmaTheta"][i + 3, i + 3]) - 1) * 100
        k[4, 0], k[4, 1] = LogMap.f_logistic_mapping(post["muTheta"][7], post["SigmaTheta"][7, 7], 2)
        # Correlation missing

    posterior.update({"k": k})

    out.update({"posterior": posterior})

    # Rap for possible next meal
    Rap = get_ModelOut(post, out_TB, ti + t[0, -1], GA_fun)[1]
    Model_Output = {"Rap": Rap - u[2, :]}

    # Model Output
    X, Ra, SigX, SigRa = get_ModelOut(post, out_TB, ti, GA_fun)

    Model_Output.update({"t": ti,
                         "G": X[[0], :],
                         "X": X[[1], :],
                         "SigG": SigX[[0], :],
                         "SigX": SigX[[1], :],
                         "Ra": Ra,
                         "SigRa": SigRa})

    out.update({"Model_Output": Model_Output})

    # Model performance
    G_pred = out_TB["suffStat"]["gx"]
    G_dat = data["y"]
    wres = (G_dat - G_pred) / (const["measCV"]/100*G_dat)
    RMSE = (G_dat - G_pred) ** 2
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


def get_ModelOut(post, out_TB, t, GA_fun):

    X = out_TB["ModelOut"]["muX"]
    n = len(out_TB["ModelOut"]["SigmaX"])
    n_theta = out_TB["options"]["dim"]["n_theta"]

    SigX = np.zeros((2, n))
    Ra = np.zeros((1, n))
    SigRa = np.zeros((1, n))

    for i in range(0, n):
        SigX[0, i] = np.sqrt(out_TB["ModelOut"]["SigmaX"][i][0, 0])
        SigX[1, i] = np.sqrt(out_TB["ModelOut"]["SigmaX"][i][1, 1])

        if GA_fun == "RaPL":
            Ra[0, i] = f_RaPL.f_Ra(t[0, i], out_TB["options"]["inF"]["tb"],
                                  post["muTheta"], out_TB["options"]["inF"]["A"],
                                  out_TB["options"]["inF"]["alpha"]) + out_TB["data"]["u"][2, i]
            dFdTh = f_RaPL.f_model(X[:, [i]], post["muTheta"], out_TB["data"]["u"][:, [i]], out_TB["options"]["inF"])[2]

        if GA_fun == "RaLN":
            Ra[0, i] = f_RaLN.f_Ra(t[0, i], post["muTheta"], out_TB["options"]["inF"]["A"])[0] + out_TB["data"]["u"][2, i]
            dFdTh = f_RaLN.f_model(X[:, [i]], post["muTheta"], out_TB["data"]["u"][:, [i]], out_TB["options"]["inF"])[2]

        dFdTh = dFdTh[0, 2:n_theta - 1] * out_TB["options"]["inF"]["V"]
        SigRa[0, i] = np.sqrt(dFdTh @ post["SigmaTheta"][2:n_theta - 1, 2:n_theta - 1] @ dFdTh.T)

    return X, Ra, SigX, SigRa


def create_ResultsFigure(out):

    opt = out["options"]
    post = out["posterior"]
    pr = out["priors"]
    data = out["data"]
    mo = out["Model_Output"]
    GA_fun = out["options"]["GA_fun"]

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 6))
    plt.subplots_adjust(left=0.08, bottom=0.05, right=0.95, top=0.95, wspace=0.3, hspace=0.3)

    # Summary
    ax[0, 0].set_title('SUMMARY')
    ax[0, 0].set_xticks([])
    ax[0, 0].set_yticks([])

    strgs = ["Results of the OMM using " + opt["GA_fun"]]
    strgs.append("")
    strgs.append("Elapsed Time " + str(out["Toolbox"]["out"]["dt"]) + " seconds")
    strgs.append("")
    strgs.append("Model Fit criteria")
    strgs.append("   - R2: " + str(np.round(out["Performance"]["R2"], 3)))
    strgs.append("   - RMSE: " + str(np.round(out["Performance"]["RMSE"], 2)) + " mmol/L")
    strgs.append("")
    strgs.append("Parameters [median +/- CV [%]]")
    strgs.append("   - p1 [1E-3 1/min]: " + str(np.round(post["p1"][0]*1E3, 2)) + " +/- "
                 + str(np.round(post["p1"][1], 1)))
    strgs.append("   - p2 [1E-3 1/min]: " + str(np.round(post["p2"][0]*1E3, 2)) + " +/- "
                 + str(np.round(post["p2"][1], 1)))
    strgs.append("   - SI [1E-4 1/min per IU]: " + str(np.round(post["SI"][0]*1E4, 2)) + " +/- "
                 + str(np.round(post["SI"][1], 1)))

    for i in range(0, len(strgs)):
        ax[0, 0].text(0.05, 0.9-i/13, strgs[i], color='k')

    # Data vs Model
    ax[0, 1].set_title('Data vs Model Prediction')
    ax[0, 1].set_ylabel("Glucose [mmol/L]")
    ax[0, 1].set_xlabel("Time [min]")

    # Data
    colors = pl.cm.Set1(np.arange(0, 2))
    ax[0, 1].plot(To1D(data["t"]), To1D(data["G"]), marker='o', markersize=3, ls='--', linewidth=0.5, color=colors[1])
    ax[0, 1].plot(To1D(data["t"]), out["const"]["Gb"]*np.ones(np.size(data["t"])), ls='--', linewidth=1, color=colors[1])
    # Model Pred
    ax[0, 1].plot(To1D(mo["t"]), To1D(mo["G"]), color=colors[1])
    ax[0, 1].fill_between(To1D(mo["t"]), To1D(mo["G"] - mo["SigG"]), To1D(mo["G"] + mo["SigG"]), alpha=0.2, color=colors[1])

    # Glucose appearance
    ax[0, 2].set_title('Glucose Appearance')
    ax[0, 2].set_ylabel("GA [mmol/kg/min]")
    ax[0, 2].set_xlabel("Time [min]")

    ax[0, 2].plot(To1D(mo["t"]), To1D(mo["Ra"]), color=colors[1])
    ax[0, 2].fill_between(To1D(mo["t"]), To1D(mo["Ra"] - mo["SigRa"]), To1D(mo["Ra"] + mo["SigRa"]), alpha=0.2, color=colors[1])

    # Weighted Residuals
    ax[1, 0].set_title('Weighted Residuals')
    ax[1, 0].set_ylabel("WRES")
    ax[1, 0].set_xlabel("Time [min]")

    ax[1, 0].plot(To1D(data["t"]), 0 * np.ones(np.size(data["t"])), ls='-', linewidth=1,
                  color='k')
    ax[1, 0].plot(To1D(data["t"]), -1 * np.ones(np.size(data["t"])), ls='--', linewidth=1,
                  color='k')
    ax[1, 0].plot(To1D(data["t"]), 1 * np.ones(np.size(data["t"])), ls='--', linewidth=1,
                  color='k')
    ax[1, 0].plot(To1D(data["t"]), To1D(out["Performance"]["wres"]), marker='o', markersize=3, ls='--', linewidth=1,
                  color=colors[1])

    # System Parameters
    ax[1, 1].set_title('System Parameters')
    ax[1, 1].set_xticks([1, 2, 3])
    ax[1, 1].set_xticklabels(['p1', 'p2', 'SI'])

    ax[1, 1].errorbar(0.9, float(pr["p1"][0] * 1E3), yerr=getLNbounds(pr["p1"], 1E3), marker='o', color=colors[0])
    ax[1, 1].errorbar(1.1, float(post["p1"][0]*1E3), yerr=getLNbounds(post["p1"], 1E3), marker='o', color=colors[1])

    ax[1, 1].errorbar(1.9, float(pr["p2"][0] * 1E3), yerr=getLNbounds(pr["p2"], 1E3), marker='o', color=colors[0])
    ax[1, 1].errorbar(2.1, float(post["p2"][0]*1E3), yerr=getLNbounds(post["p2"], 1E3), marker='o', color=colors[1])

    ax[1, 1].errorbar(2.9, float(pr["SI"][0] * 1E4), yerr=getLNbounds(pr["SI"], 1E4), marker='o', color=colors[0])
    ax[1, 1].errorbar(3.1, float(post["SI"][0]*1E4), yerr=getLNbounds(post["SI"], 1E4), marker='o', color=colors[1])

    # Input Parameters
    ni = np.shape(post["k"])[0]
    ax[1, 2].set_title('Input Parameters')
    ax[1, 2].set_xticks(range(1, ni+1))

    if GA_fun == "RaPL":
        for i in range(1, ni+1):
            ax[1, 2].errorbar(i-0.1, float(pr["k"][i-1, 0] * 1), yerr=getLNbounds(pr["k"][i-1, :], 1), marker='o', color=colors[0])
            ax[1, 2].errorbar(i+0.1, float(post["k"][i-1, 0] * 1), yerr=getLNbounds(post["k"][i-1, :], 1), marker='o',
                              color=colors[1])
            ax[1, 2].set_xticklabels(['k1', 'k2', 'k3', 'k4', 'k5', 'k7'])

    if GA_fun == "RaLN":
        sc = [0.1, 1, 0.1, 1]
        for i in range(1, ni):
            ax[1, 2].errorbar(i-0.1, float(pr["k"][i-1, 0] * sc[i-1]), yerr=getLNbounds(pr["k"][i-1, :], sc[i-1]), marker='o',
                              color=colors[0])
            ax[1, 2].errorbar(i+0.1, float(post["k"][i-1, 0] * sc[i-1]), yerr=getLNbounds(post["k"][i-1, :], sc[i-1]), marker='o',
                              color=colors[1])

        ax[1, 2].errorbar(5 - 0.1, float(pr["k"][4, 0] * 1), yerr=getLogisticbounds(pr["k"][4, :]), marker='o',
                          color=colors[0])
        ax[1, 2].errorbar(5 + 0.1, float(post["k"][4, 0] * 1), yerr=getLogisticbounds(post["k"][4, :]), marker='o',
                          color=colors[1])
        ax[1, 2].set_xticklabels(['T1', 'W1', 'T2', 'W2', 'Rh'])


    plt.show()


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