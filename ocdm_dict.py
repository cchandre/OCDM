###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp

Method = 'trajectories'

dimension = 2

E0 = 2e-2
Omega = lambda t: 7e-7 * t
envelope = 'sinus'
te = [10, 50, 10]

Energy0 = -0.07

r = [2.5, 10]
Ntraj = 10
type_traj = ['all', 'spherical', 'rotated']
dpi = 3000

SaveData = False
PlotResults = True
Parallelization = (False, 50)

darkmode = True

###################################################################################################
##                             DO NOT EDIT BELOW                                                 ##
###################################################################################################
if Method == 'poincaré':
    dimension = 2
    envelope = 'const'
    te = [0, sum(te), 0]
    type_traj = ['all', 'spherical', 'rotated']
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
        'Tol': [1e-9, 1e-7],
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'contour_levels': 50,
        'darkmode': darkmode})
###################################################################################################
