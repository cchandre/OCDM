###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp

Method = 'poincaré'

dimension = 2

E0 = 5e-2
Omega = lambda t: 2e-5
envelope = 'sinus'
te = [15, 100, 15]

Energy0 = -0.001

r = [2.5, 15]
Ntraj = 50
type_traj = ['dissociated', 'spherical', 'rotated']
dpi = 3000

SaveData = False
PlotResults = True
Parallelization = (False, 50)

darkmode = True

###################################################################################################
##                             DO NOT EDIT BELOW                                                 ##
###################################################################################################
if Method == 'poincaré':
    envelope = 'const'
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
        'type_traj': type_traj,
        'ode_solver': 'RK45',
        'Tol': [1e-8, 1e-5],
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'contour_levels': 50,
        'darkmode': darkmode})
###################################################################################################
