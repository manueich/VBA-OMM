import numpy as np
import Functions.VBA_Python.VBA as VBA
import Functions.VBA_Python.VBA_NoiseDistConv as NoiseDistConv
import Functions.logistic_mapping as LogMap
import Functions.f_OMM_RaPL as f_RaPL
import Functions.f_OMM_RaLN as f_RaLN

def VBA_OMM_G(dat, priors, const, opt):

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

    # Check const
    try:
        A = const["A"]
        V = const["V"]
        dt = const["dt"]
        if not const["Rap"]:
            Rap = np.zeros((1, int(t[0, -1]/dt+1)))
        else:
            Rap = np.reshape(const["Rap"],(1,const["Rap"]))
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

    data = {"y": G[:, 1:np.size(G)],
            "t": t[:, 1:np.size(G)],
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
    for i in range(1, np.size(t) - 1):
        iQy.append(np.eye(1)*1/(G[0, i+1]/np.mean(Gi)))

    a, b = NoiseDistConv.ToGamma(np.mean(Gi)*measCV/100, np.mean(Gi)*measCV/100*0.1)

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
               "verbose": True,
               "Display": False,
               "updateHP": False,
               "TolFun": 1E-4,
               "GnTolFun": 1e-4,
               "MaxIter": 100}

    # -----------------------------------------
    # ---- INVERSION --------

    if displayWin:
        print("Model Inversion ...")

    #try:
    post, out_TB = VBA.main(data, ti, pr, options)
    #except:
    print("Inverison Failed")
    #return []

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

    return out
