import numpy as np
from Functions.VBA_Python import VBA_UpdateHP as VBA_HP, VBA_GaussNewton as VBA_GN, VBA_Initialise as VBA_init, \
    VBA_plotting as plotting, VBA_basics as base
from timeit import default_timer as timer
import matplotlib.pyplot as plt


def main(data, t, priors, options):

    # Check inputs Get initial estimates with priors
    posterior, priors, options, suffStat = VBA_init.Initialize(data, t, priors.copy(), options.copy())

    # Plot data and priors
    if options["Display"]:
        ax = plotting.plot_data(data)
        ax = plotting.plot_model(ax, suffStat, posterior, priors, data, options)

    # -----------------------------------------------------------
    # Main iteration loop maximising Free Energy
    stop = False
    it = 0
    start = timer()
    if options["verbose"]:
        print("Maximising Free Energy ...")
        print('Prior Free Energy         F=', np.round(float(suffStat["F"][-1]), 1))

    while not stop:
        it = it + 1
        F0 = suffStat["F"][-1]

        # Update Evolution parameters
        posterior, suffStat = VBA_GN.GaussNewton(data, t, posterior.copy(), priors, suffStat.copy(), options)

        # Update Noise parameters
        if options["updateHP"]:
            posterior, suffStat = VBA_HP.UpdateHP(data, t, posterior.copy(), priors, suffStat.copy(), options)

        # Display progress
        if options["verbose"]:
            dF = suffStat["F"][-1] - F0
            print('VB iteration #', it, '         F=', np.round(float(suffStat["F"][-1]), 1), '         ... dF=', np.round(float(dF), 3))

        # Plot current State
        if options["Display"]:
            ax = plotting.plot_model(ax, suffStat, posterior, priors, data, options)

        # Check Convergence Criterion
        F = suffStat["F"]
        dF = F[-1] - F[-2]

        if it == options["MaxIter"] or np.abs(dF) <= options["TolFun"]:
            stop = True
            if np.abs(dF) <= options["TolFun"]:
                conv = 1
            else:
                conv = 0
    end = timer()
    print('VB inversion complete (took ~', np.round(end-start, 1), 's).')

    del posterior["iQy"]
    del posterior["muP"]
    del posterior["SigmaP"]

    out = {"it": it,
           "F": F[-1],
           "priors": priors,
           "ModelOut": suffStat["model_out"],
           "suffStat": suffStat}

    # # Keep Figure Window displayed
    # if options["Display"]:
    #     plt.show()

    return posterior, out


def simulate(td, t, u, priors, options, plot_bool):

    # Check inputs
    options = VBA_init.check_options(options.copy())
    priors = VBA_init.check_priors(options, priors.copy())
    data, options, priors = check_data(td, t, u, priors.copy(), options.copy())
    options, priors = check_model(options.copy(), priors.copy(), data)

    # Get initial estimates with priors
    try:
        y, muX, SigmaX = base.solveODE(t, priors["muP"], priors["SigmaP"], data["u"], options)[0:3]
    except:
        raise Exception("The model produces error (see above)")

    if base.isWeird(y):
        raise Exception("Could not simulate model: model generates NaN or Inf!'")

    yd = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))
    epsilon = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))

    for i in range(0, td.size):
        idx = np.where(np.round(t, 8) == np.round(td[0, i], 8))

        sigma = priors["a"] / priors["b"]
        iQyt = priors["iQy"][i]

        C = base.Invert_M(iQyt * sigma)
        epsilon[:, [i]] = np.reshape(np.random.multivariate_normal(np.zeros(options["dim"]["nY"]), C),
                                     (options["dim"]["nY"], 1))

        yd[:, [i]] = y[:, idx[1]] + epsilon[:, [i]]

    if plot_bool:
        plotting.plot_sim(td, yd, t, muX)

    return yd


def compare_to_sim(posterior, prios_sim, options):

    plotting.compare_to_sim(posterior, prios_sim, options)
    plt.show()

    return


def check_data(ty, t, u, priors, options):
    dim = options["dim"]

    # Check ty
    if np.shape(ty)[0] != 1:
        raise Exception("ty must be a 1 by n array")
    else:
        data = ({"ty": ty})
        dim.update({"nD": np.shape(ty)[1]})

    # Check u
    if isinstance(u, list):
        data.update({"u": np.zeros((1, t.size))})
        dim.update({"nu": 0})
    elif np.shape(u)[1] != np.shape(t)[1]:
        raise Exception("Inputs in u must be specified on the ODE integration time step t")
    elif base.isWeird(u):
        raise Exception("The data in u contains NaNs or Infs")
    else:
        dim.update({"nu": np.shape(u)[0]})
        data.update({"u": u})

    # Check integration time Grid
    if t[0, 0] != 0:
        raise Exception("The ODE integration time grid must begin with 0")
    dt = t[0, 1] - t[0, 0]
    if np.any(np.round(np.diff(t), 8) != dt):
        raise Exception("The ODE integration time grid must be equally spaced")
    if ty[0, 0] < 0 or ty[0, -1] > t[0, -1]:
        raise Exception("Data timepoints ty lie outside of integration time grid")

    options.update({"dim": dim})

    return data, options, priors


def check_model(options, priors, data):
    dim = options["dim"]
    muP = priors["muP"]
    u = data["u"]

    phi = muP[0: dim["n_phi"]]
    th = muP[dim["n_phi"]: dim["n_theta"] + dim["n_phi"]]
    x0 = muP[dim["n_theta"] + dim["n_phi"]: dim["n_theta"] + dim["n_phi"] + dim["n"]]

    # Get functions
    f_model = options["f_model"]
    f_obs = options["f_obs"]

    try:
        x, J, H = f_model(x0, th, u[:, [0]], options["inF"])
        y, dY_dX, dY_dPhi = f_obs(x, phi, u[:, [0]], options["inG"])
    except:
        raise Exception("The model produces error (see above)")

    if np.shape(x)[0] != dim["n"] or np.shape(x)[1] != 1:
        raise Exception("Model Error: Dimensions of x must be n by 1")

    dim.update({"nY": np.shape(y)[0]})
    if np.shape(y)[1] != 1:
        raise Exception("Model Error: Dimensions of y must be nY by 1")

    options.update({"dim": dim})

    # Check iQy now that we know nY
    if "iQy" in priors:
        if len(priors["iQy"]) != options["dim"]["nD"]:
            raise Exception("The size of iQy must match the given data")
        else:
            iQy = [0] * len(priors["iQy"])
            for i in range(0, len(priors["iQy"])):
                if np.shape(priors["iQy"][i])[0] != options["dim"]["nY"] or np.shape(priors["iQy"][i])[1] != \
                        options["dim"]["nY"]:
                    raise Exception("Inconsistent dimension in iQy")
                else:
                    diQ = np.diag(priors["iQy"][i])
                    iQy[i] = np.diag(diQ) @ priors["iQy"][i] @ np.diag(diQ)
    else:
        iQy = [np.eye(options["dim"]["nY"])]
        for i in range(0, options["dim"]["nD"] - 1):
            iQy.append(np.eye(options["dim"]["nY"]))

    priors.update({"iQy": iQy})

    return options, priors


