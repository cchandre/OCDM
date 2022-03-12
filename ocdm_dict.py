###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp

Method = 'trajectories'
r = [2.5, 15]

dimension = 2

E0 = 5e-1
Omega = lambda t: 2e-5
envelope = 'sinus'
te = [15, 70, 15]

Energy0 = -0.001

Ntraj = 2
Tol = 1e-5
plot_traj = 'all'

SaveData = False
PlotResults = False
Parallelization = (False, 3)

dpi = 1000
darkmode = True

###################################################################################################
##                             DO NOT EDIT BELOW                                                 ##
###################################################################################################
dict_list = [{'Method': Method} for _ in xp.atleast_1d(E0)]
for _, E in enumerate(xp.atleast_1d(E0)):
    dict_list[_].update({
        'r': r,
        'dim': dimension,
        'E0': E,
        'envelope': envelope,
        'te': xp.asarray(te) / 2.42e-5,
        'Energy0': Energy0,
        'Ntraj': Ntraj,
        'Tol': Tol,
		'type_traj': type_traj,
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'dpi': dpi,
        'darkmode': darkmode})
###################################################################################################
