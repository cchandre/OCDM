###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp

Method = 'dissociation'

dimension = 2

E0 = xp.linspace(0, 1e-1, 150)
Omega = lambda t: 2e-5
envelope = 'sinus'
te = [15, 50, 15]

Energy0 = -0.001

r = [2.5, 15]
Ntraj = 500
plot_traj = ['dissociated', 'spherical']
dpi = 300

SaveData = False
PlotResults = False
Parallelization = (True, 3)

darkmode = True

###################################################################################################
##                             DO NOT EDIT BELOW                                                 ##
###################################################################################################
dict_list = [{'Method': Method} for _ in xp.atleast_1d(E0)]
for dict, E in zip(dict_list, xp.atleast_1d(E0)):
    dict.update({
        'dim': dimension,
        'E0': E,
        'envelope': envelope,
        'te': xp.asarray(te) / 2.42e-5,
        'Energy0': Energy0,
        'r': r,
        'dpi': dpi,
        'Ntraj': Ntraj,
        'plot_traj': plot_traj,
        'ode_solver': 'RK45',
        'Tol': [1e-6, 1e-3],
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'contour_levels': 50,
        'darkmode': darkmode})
###################################################################################################
