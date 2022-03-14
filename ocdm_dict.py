###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp

Method = 'dissociation'
r = [2.5, 15]

dimension = 2

E0 = xp.linspace(0, 1e-1, 150)
Omega = lambda t: 2e-5
envelope = 'sinus'
te = [15, 50, 15]

Energy0 = -0.001

Ntraj = 500
Tol = [1e-5, 1e-3]
plot_traj = ['dissociated', 'spherical']

SaveData = False
PlotResults = False
Parallelization = (True, 3)

dpi = 300
darkmode = True

###################################################################################################
##                             DO NOT EDIT BELOW                                                 ##
###################################################################################################
dict_list = [{'Method': Method} for _ in xp.atleast_1d(E0)]
for dict, E in zip(dict_list, xp.atleast_1d(E0)):
    dict.update({
        'r': r,
        'dim': dimension,
        'E0': E,
        'envelope': envelope,
        'te': xp.asarray(te) / 2.42e-5,
        'Energy0': Energy0,
        'Ntraj': Ntraj,
        'Tol': Tol,
		'plot_traj': plot_traj,
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'dpi': dpi,
        'contour_levels': 50,
        'darkmode': darkmode})
###################################################################################################
