import numpy as np
from timeit import default_timer as timer
from . import VBA_Initialize
from . import VBA_plotting
from . import VBA_GaussNewton
from . import VBA_UpdateHP
from . import VBA_basics
from . import VBA_getDiagnostics
import matplotlib.pyplot as plt


def main(data, t, priors, options):

    # Check inputs Get initial estimates with priors
    posterior, priors, options, suffStat = VBA_Initialize.Initialize(data, t, priors.copy(), options.copy())

    # Plot data and priors
    if options["Display"]:
        ax = VBA_plotting.plot_data(data)
        ax = VBA_plotting.plot_model(ax, suffStat, posterior, priors, data, options)

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
        posterior, suffStat = VBA_GaussNewton.GaussNewton(data, t, posterior.copy(), priors, suffStat.copy(), options)

        # Update Noise parameters
        if options["updateHP"]:
            posterior, suffStat = VBA_UpdateHP.UpdateHP(data, t, posterior.copy(), priors, suffStat.copy(), options)

        # Display progress
        if options["verbose"]:
            dF = suffStat["F"][-1] - F0
            print('VB iteration #', it, '         F=', np.round(float(suffStat["F"][-1]), 1), '         ... dF=', np.round(float(dF), 3))

        # Plot current State
        if options["Display"]:
            ax = VBA_plotting.plot_model(ax, suffStat, posterior, priors, data, options)

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
    if options["verbose"]:
        print('VB inversion complete (took ~', np.round(end-start, 1), 's).')

    del posterior["iQy"]
    del posterior["muP"]
    del posterior["SigmaP"]

    out = {"it": it,
           "F": F[-1],
           "priors": priors,
           "ModelOut": suffStat["model_out"],
           "suffStat": suffStat,
           "options": options,
           "data": data,
           "dt": np.round(end-start, 1)}

    fit = VBA_getDiagnostics.get_fit(data, priors, posterior, options, suffStat)

    out.update({"fit": fit})

    # Keep Figure Window displayed
    if options["Display"]:
        plt.show()

    return posterior, out


