###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp

Method = 'display_potentials'
r = [2.5, 15]

dimension = 2

E0 = 5e-3
Omega = lambda t: 0
envelope = 'const'
te = [15, 70, 15]

Energy0 = -0.01

Ntraj = 100
Tol = 1e-3

SaveData = False
PlotResults = True

dpi = 300
darkmode = True

###################################################################################################
##                             DO NOT EDIT BELOW                                                 ##
###################################################################################################
dict = {
        'Method': Method,
        'r': r,
        'dim': dimension,
        'E0': E0,
        'Omega': Omega,
        'envelope': envelope,
        'te': xp.asarray(te) / 2.42e-5,
        'Energy0': Energy0,
        'Ntraj': Ntraj,
        'Tol': Tol,
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'dpi': dpi,
        'darkmode': darkmode}
###################################################################################################
