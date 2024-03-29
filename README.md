# VBA-OMM
Identification of the Oral Minimal Model of glucose dynamcis from non-fasting conditions using Variational Bayesian Analysis.

**UPDATE: The Python implementation now also contains the OMM of C-peptide dynamics to infer beta-cell responsivity. See the respective *__init__.py* (in folder VBA-OMM/Python/VBA_OMM/) script for details.**

# MATLAB
To use the MATLAB version of the toolbox simiply download all files and subfolders in folder MATLAB and add them to your path. The main function is called VBA_OMM_G.m and contains a detailed description on how to use it. To see how the function is used and its inputs are defined, you can run the respective script in the \demo folder. A detailed description of the underlying methodology is found in the open-access publication below. The core of the VBA approach in MATLAB is copied from the following wesite http://mbb-team.github.io/VBA-toolbox/.

# Python 
To use the Python version of the toolbox is in the folder called Python. It contains a package called "VBA_OMM" with a detailed description in the __init__.py file. The package contains an implementation of the VBA method from my other repository (https://github.com/manueich/VBA-python). To see how the package is used look at the script demo_OMM.py. A detailed description of the underlying methodology is found in the open-access publication below.

# Structural Identifiability
A structural indetifiability analysis of the OMM using the Taylor series approach implemented in Mathematica

# About
This software has been developed in the context of a PhD project under the supervison of Natasha Khovanova from the School of Engineering at the University of Warwick and John Hattersley from the University Hospitals Coventry and Warwickshire.

If you have any questions you can contact me under manuel.eichenlaub@gmail.com

# License
This software is distributed under a GNU open-source licence. You are free to download and modify it, but we would be most grateful, should you acknowledge the authors by citing the following reference:

M. Eichenlaub, J. Hattersley, M. Gannon, F. Nuttall, N. Khovanova (2021): "Bayesian parameter estimation in the oral minimal model of glucose dynamics from non-fasting conditions using a new function of glucose appearance". In: Computer Methods and Programs in Biomedicine, 200 (https://doi.org/10.1016/j.cmpb.2020.105911)

Additionally, it would be appreciated if you would acknowledge the authors of the underlying variational Bayesian appraoch by citing the following reference:

J. Daunizeau, V. Adam, L. Rigoux (2014), VBA: a probabilistic treatment of nonlinear models for neurobiological and behavioural data. PLoS Comp Biol 10(1): e1003441.

M. Eichenlaub 05/01/2021
