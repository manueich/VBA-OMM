import numpy as np
import scipy.linalg as la
from . import VBA_basics as base


def Initialize(data, t, priors, options):

    # Check inputs
    options = check_options(options.copy())
    priors = check_priors(options, priors.copy())
    options, priors = check_data(options.copy(), priors.copy(), data, t)
    check_model(options, priors, data)

    posterior = priors.copy()

    # Get initial estimates with priors
    try:
        y, muX, SigmaX, dXdTh, dXdX0, dYdPhi, dYdTh, dYdX0, dG_dP = base.solveODE(t, priors["muP"], priors["SigmaP"],
                                                                                data["u"], options)
    except:
        raise Exception("The model produces error (see above)")


    if base.isWeird(y):
        raise Exception("Could not initialize VB scheme: model generates NaN or Inf!'")

    yd = data["y"]
    td = data["t"]

    dy = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))
    vy = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))
    gx = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))
    dy2 = 0
    logL = 0
    div = False

    # Loop over time series
    for i in range(0, td.size):
        idx = np.where(np.round(t, 8) == np.round(td[0, i], 8))

        sigmaHat = priors["a"] / priors["b"]
        iQyt = priors["iQy"][i]

        # Get intermediate values
        gx[:, [i]] = y[:, idx[1]]
        _, _, logL_t, dy[:, [i]], dy2_t, vy[:, [i]] = get_dL(y[:, idx[1]], yd[:, [i]], dG_dP[int(idx[1])], sigmaHat, iQyt)

        # Aggregate
        dy2 = dy2 + dy2_t
        logL = logL + logL_t

        V = dG_dP[int(idx[1])].T @ priors["SigmaP"] @ dG_dP[int(idx[1])]
        # print(V)
        vy[:, [i]] = vy[:, [i]] + (np.diag(V)).reshape((options["dim"]["nY"], 1))

        if base.isWeird(dy2) or base.isWeird(dG_dP[int(idx[1])]):
            div = True
            break

    # Collect results
    model_out = {"t": t,
                 "y": y,
                 "muX": muX,
                 "SigmaX": SigmaX,
                 "dXdTh": dXdTh,
                 "dXdX0": dXdX0,
                 "dYdPhi": dYdPhi,
                 "dYdTh": dYdTh,
                 "dYdX0": dYdX0,
                 "dG_dP": dG_dP}

    suffStat = {"gx": gx,
                "dy": dy,
                "vy": vy,
                "dy2": dy2,
                "logL": logL,
                "dP": np.zeros((options["dim"]["n_phi"] + options["dim"]["n_theta"] + options["dim"]["n"], 1)),
                "SigmaP": priors["SigmaP"],
                "model_out": model_out,
                "data": data}

    F = base.Free_Energy(posterior, priors, suffStat, options)

    suffStat.update({"F": [F]})

    return posterior, priors, options, suffStat


def check_options(options):
    if not "dim" in options:
        raise Exception("Please provide model dimensions dim")

    dim = options["dim"]

    # Check input dimensions
        # - Model Order
    if not "n" in dim:
        raise Exception("Please provide model order n in dim")
    elif dim["n"] < 0 or (dim["n"] % 1) != 0:
        raise Exception("Please provide a valid model order n in dim")
    if dim["n"] == 0:
        dim.update({"n_theta": 0})
        # - No of evolution parameters
    if not "n_theta" in dim:
        raise Exception("Please provide the number of evolution parameter n_theta in dim")
    elif dim["n_theta"] < 0 or (dim["n_theta"] % 1) != 0:
        raise Exception("Please provide a valid number of evolution parameters n_theta in dim")
        # - No of observation parameters
    if not "n_phi" in dim:
        raise Exception("Please provide the number of observation parameters in dim")
    elif dim["n_phi"] < 0 or (dim["n_phi"] % 1) != 0:
        raise Exception("Please provide a valid number of observation parameters n_phi in dim")

    if not "f_obs" in options:
        raise Exception("Please provide observation function")
    if not "f_model" in options:
        raise Exception("Please provide model function")

    options.update({"dim": dim})

    # Check other options
    if not "GnMaxIter" in options or options["GnMaxIter"] <= 0:
        options.update({"GnMaxIter": 32})
    if not "GnTolFun" in options or options["GnTolFun"] <= 0:
        options.update({"GnTolFun": 1e-5})
    if not "MaxIter" in options or options["MaxIter"] <= 0:
        options.update({"MaxIter": 32})
    if not "TolFun" in options or options["TolFun"] <= 0:
        options.update({"TolFun": 1e-5})
    if not "updateHP" in options:
        options.update({"updateHP": True})
    if not "verbose" in options:
        options.update({"verbose": True})
    if not "Display" in options:
        options.update({"Display": True})
    if not "ODESolver" in options:
        options.update({"ODESolver": 'Euler'})
    elif options["ODESolver"] != 'Euler':
        if options["ODESolver"] != 'RK':
            raise Exception("Please specify either 'RK' or 'Euler' as ODESolver in options")
    if not "inF" in options:
        options.update({"inF": []})
    if not "inG" in options:
        options.update({"inG": []})

    return options


def check_priors(options, priors):
    # Check priors for initial conditions
    if options["dim"]["n"] > 0:
        if not "muX0" in priors or not "SigmaX0" in priors:
            raise Exception("Please specify priors for X0")
        elif options["dim"]["n"] != np.size(priors["muX0"]):
            raise Exception("Dimension of priors for muX0 does not match model order")
        elif not np.allclose(options["dim"]["n"] * np.ones((1, 2)), np.shape(priors["SigmaX0"]), rtol=1e-12,
                             atol=1e-12):
            raise Exception("Dimension of priors for SigmaX0 does not match model order")
        if np.any(np.diag(priors["SigmaX0"]) == 0):
            idx = np.where(np.diag(priors["SigmaX0"]) == 0)
            sig = priors["SigmaX0"]
            for i in range(0, np.size(idx)):
                sig[idx, idx] = 1E-12
            priors.update({"SigmaX0": sig})
    else:
        priors.update({"muX0": np.array([[]]).T, "SigmaX0": np.array([[]]).T})

    # Check priors for evolution Parameters
    if options["dim"]["n_theta"] > 0:
        if not "muTheta" in priors or not "SigmaTheta" in priors:
            raise Exception("Please specify priors for Theta")
        elif options["dim"]["n_theta"] != np.size(priors["muTheta"]):
            raise Exception("Dimension of priors for muTheta does not match n_theta in dim")
        elif not np.allclose(options["dim"]["n_theta"] * np.ones((1, 2)), np.shape(priors["SigmaTheta"]), rtol=1e-12,
                             atol=1e-12):
            raise Exception("Dimension of priors for SigmaTheta does not match n_theta in dim")
        if np.any(np.diag(priors["SigmaTheta"]) == 0):
            idx = np.where(np.diag(priors["SigmaTheta"]) == 0)
            sig = priors["SigmaTheta"]
            for i in range(0, np.size(idx)):
                sig[idx, idx] = 1e-12
            priors.update({"SigmaTheta": sig})
    else:
        priors.update({"muTheta": np.array([[]]).T, "SigmaTheta": np.array([[]]).T})

    # Check priors for observation Parameters
    if options["dim"]["n_phi"] > 0:
        if not "muPhi" in priors or not "SigmaPhi" in priors:
            raise Exception("Please specify priors for Phi")
        elif options["dim"]["n_phi"] != np.size(priors["muPhi"]):
            raise Exception("Dimension of priors for muPhi does not match n_phi in dim")
        elif not np.allclose(options["dim"]["n_phi"] * np.ones((1, 2)), np.shape(priors["SigmaPhi"]), rtol=1e-12,
                             atol=1e-12):
            raise Exception("Dimension of priors for SigmaPhi does not match n_phi in dim")
        if np.any(np.diag(priors["SigmaPhi"]) == 0):
            idx = np.where(np.diag(priors["SigmaPhi"]) == 0)
            sig = priors["SigmaPhi"]
            for i in range(0, np.size(idx)):
                sig[idx, idx] = 1e-12
            priors.update({"SigmaPhi": sig})
    else:
        priors.update({"muPhi": np.array([[]]).T, "SigmaPhi": np.array([[]]).T})

    # Concatenation
    if options["dim"]["n_phi"] > 0:
        SigmaP = priors["SigmaPhi"]
        muP = priors["muPhi"]
        if options["dim"]["n"] > 0:
            if options["dim"]["n_theta"] > 0:
                SigmaP = la.block_diag(SigmaP, priors["SigmaTheta"])
                muP = np.concatenate((muP, priors["muTheta"]), 0)
                SigmaP = la.block_diag(SigmaP, priors["SigmaX0"])
                muP = np.concatenate((muP, priors["muX0"]), 0)
            else:
                SigmaP = la.block_diag(SigmaP, priors["SigmaX0"])
                muP = np.concatenate((muP, priors["muX0"]), 0)
    else:
        if options["dim"]["n_theta"] > 0:
            SigmaP = priors["SigmaTheta"]
            muP = priors["muTheta"]
            SigmaP = la.block_diag(SigmaP, priors["SigmaX0"])
            muP = np.concatenate((muP, priors["muX0"]), 0)
        else:
            SigmaP = priors["SigmaX0"]
            muP = priors["muX0"]

    priors.update({"muP": muP, "SigmaP": SigmaP})

    #  Check Noise priors
    if "a" in priors:
        if priors["a"] <= 0:
            priors["a"] = 1
            print('Warning: Prior value of a was negative. Default value of 1 was set.')
    else:
        priors["a"] = 1
        print('Warning: Prior value of a was missing. Default value of 1 was set.')

    if "b" in priors:
        if priors["b"] <= 0:
            priors["b"] = 1
            print('Warning: Prior value of b was negative. Default value of 1 was set.')
    else:
        priors["b"] = 1
        print('Warning: Prior value of b was missing. Default value of 1 was set.')


    return priors


def check_data(options, priors, data, t):

    if not "y" in data:
        raise Exception("Please provide the data in y")
    else:
        y = data["y"]
    if not "t" in data:
        raise Exception("Please provide the time for data y")
    else:
        ty = data["t"]

    dim = options["dim"]

    # Check Data
    if np.shape(ty)[0] != 1:
        raise Exception("ty must be a 1 by nD array")
    if np.shape(ty)[1] != np.shape(y)[1]:
        raise Exception("The length of ty must match the elements in y")
    if base.isWeird(y):
        raise Exception("The data in y contains NaNs or Infs")

    # Check Input
    if not "u" in data:
        # data.update({"u": np.zeros((1, t.size))})
        dim.update({"nu": 0})
    elif np.shape(data["u"])[1] != np.shape(t)[1]:
        raise Exception("Inputs in u must be specified on the ODE integration time step t")
    elif base.isWeird(data["u"]):
        raise Exception("The data in u contains NaNs or Infs")
    else:
        dim.update({"nu": np.shape(data["u"])[0]})

    # Check integration time Grid
    if t[0, 0] != 0:
        raise Exception("The ODE integration time grid must begin with 0")
    dt = t[0, 1] - t[0, 0]
    if np.any(np.round(np.diff(t), 8) != dt):
        raise Exception("The ODE integration time grid must match dt in inF")
    if ty[0, 0] < 0 or ty[0, -1] > t[0, -1]:
        raise Exception("Data timepoints ty lie outside of integration time grid")

    dim.update({"nD": np.shape(y)[1]})
    dim.update({"nY": np.shape(y)[0]})
    options.update({"dim": dim})

    # Check iQy
    if "iQy" in priors:
        if len(priors["iQy"]) != options["dim"]["nD"]:
            raise Exception("The size of iQy must match the given data")
        else:
            for i in range(0, len(priors["iQy"])):
                if np.shape(priors["iQy"][i])[0] != options["dim"]["nY"] or np.shape(priors["iQy"][i])[1] != options["dim"]["nY"]:
                    raise Exception("Inconsistent dimension in iQy")
    else:
        iQy = [np.eye(options["dim"]["nY"])]
        for i in range(0, options["dim"]["nD"]-1):
            iQy.append(np.eye(options["dim"]["nY"]))

    return options, priors


def check_model(options, priors, data):
    # Check the model functions

    dim = options["dim"]
    muP = priors["muP"]
    u = data["u"]

    phi = muP[0: dim["n_phi"]]
    th = muP[dim["n_phi"]: dim["n_theta"]+dim["n_phi"]]
    x0 = muP[dim["n_theta"]+dim["n_phi"]: dim["n_theta"]+dim["n_phi"]+dim["n"]]

    # Get functions
    f_model = options["f_model"]
    f_obs = options["f_obs"]

    try:
        dx, J, H = f_model(x0, th, u[:, [0]], options["inF"])
        y, dY_dX, dY_dPhi = f_obs(dx, phi, u[:, [0]], options["inG"])
    except:
        raise Exception("The model produces error (see above)")

    if np.shape(dx)[0] != dim["n"] or np.shape(dx)[1] != 1:
        raise Exception("Model Error: Dimensions of x must be n by 1")
    if np.shape(J)[0] != dim["n"] or np.shape(J)[1] != dim["n"]:
        raise Exception("Model Error: Dimensions of J must be n by n")
    if np.shape(H)[0] != dim["n"] or np.shape(H)[1] != dim["n_theta"]:
        raise Exception("Model Error: Dimensions of H must be n by n_theta")

    if np.shape(y)[0] != dim["nY"] or np.shape(y)[1] != 1:
        raise Exception("Model Error: Dimensions of y must be nY by 1")
    if not dY_dX.size != 0:
        if np.shape(dY_dX)[0] != dim["nY"] or np.shape(dY_dX)[1] != dim["n"]:
            raise Exception("Model Error: Dimensions of dY_dX must be 1 by n")
    if not dY_dPhi.size != 0:
        if np.shape(dY_dPhi)[0] != dim["nY"] or np.shape(dY_dPhi)[1] != dim["n_phi"]:
            raise Exception("Model Error: Dimensions of dY_dPhi must be 1 by n_phi")


def get_dL(y, yd, dG_dP, sigmaHat, iQyt):
    # Compute useful intermediate values describing the misfit between a model
    # prediction and an observation.
    #
    # IN:
    #   - gx: model prediction about the observation y (1st order moment)
    #   - dG_dPhi: derivative of gx wrt observation parameters
    #   - y: actual observation
    #   - sigmaHat: scaling factor of the 2nd order moment of the prediction
    #
    # OUT:
    #   - ddydP: gradient of the prediction error wrt observation parameters
    #   - d2gdx2: hessian of the prediction
    #   - logL: log-likelihood of the observation given the prediction
    #   - dy: prediction error
    #   - dy2: normalized squared deviation
    #   - vy: prediction variance

    dy = yd - y

    vy = 1 / sigmaHat * np.diag(base.Invert_M(iQyt))
    vy = vy.reshape((np.shape(iQyt)[0], 1))
    ddydP = sigmaHat * (dG_dP @ iQyt @ dy)
    d2gdx2 = sigmaHat * (dG_dP @ iQyt @ dG_dP.T)
    dy2 = dy.T @ iQyt @ dy
    logL = -0.5 * sigmaHat * dy2 + 0.5 * base.log_det(iQyt *sigmaHat) - 0.5 * np.log(2 * np.pi)

    return ddydP, d2gdx2, logL, dy, dy2, vy
