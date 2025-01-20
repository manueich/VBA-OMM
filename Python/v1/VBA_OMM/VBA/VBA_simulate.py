import numpy as np
from . import VBA_Initialize
from . import VBA_plotting
from . import VBA_basics
import matplotlib.pyplot as plt


def simulate(td, t, u, priors, options, plot_bool):

    # Check inputs
    options = VBA_Initialize.check_options(options.copy())
    priors = VBA_Initialize.check_priors(options, priors.copy())
    data, options, priors = check_data(td, t, u, priors.copy(), options.copy())
    options, priors = check_model(options.copy(), priors.copy(), data)

    # Get initial estimates with priors
    try:
        y, muX, SigmaX = VBA_basics.solveODE(t, priors["muP"], priors["SigmaP"], data["u"], options)[0:3]
    except:
        raise Exception("The model produces error (see above)")

    if VBA_basics.isWeird(y):
        raise Exception("Could not simulate model: model generates NaN or Inf!'")

    yd = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))
    epsilon = np.zeros((options["dim"]["nY"], options["dim"]["nD"]))

    for i in range(0, td.size):
        idx = np.where(np.round(t, 8) == np.round(td[0, i], 8))

        sigma = priors["a"] / priors["b"]
        iQyt = priors["iQy"][i]

        C = VBA_basics.Invert_M(iQyt * sigma)
        epsilon[:, [i]] = np.reshape(np.random.multivariate_normal(np.zeros(options["dim"]["nY"]), C),
                                     (options["dim"]["nY"], 1))

        yd[:, [i]] = y[:, idx[1]] + epsilon[:, [i]]

    if plot_bool:
        VBA_plotting.plot_sim(td, yd, t, muX)

    return yd


def compare_to_sim(posterior, prios_sim, options):

    VBA_plotting.compare_to_sim(posterior, prios_sim, options)
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
    elif VBA_basics.isWeird(u):
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
            for i in range(0, len(priors["iQy"])):
                if np.shape(priors["iQy"][i])[0] != options["dim"]["nY"] or np.shape(priors["iQy"][i])[1] != \
                        options["dim"]["nY"]:
                    raise Exception("Inconsistent dimension in iQy")
    else:
        iQy = [np.eye(options["dim"]["nY"])]
        for i in range(0, options["dim"]["nD"] - 1):
            iQy.append(np.eye(options["dim"]["nY"]))

    priors.update({"iQy": iQy})

    return options, priors


