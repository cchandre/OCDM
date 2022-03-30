###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp

Method = 'dissociation'

dimension = 2

E0 = xp.linspace(0, 0.1, 30)
Omega = lambda t: 3e-10 * t
envelope = 'sinus'
te = [5, 40, 5]

Ntraj = 100
r = [2.5, 10]
initial_J = 30
EnergyPS = []

type_traj = ['dissociated', 'spherical', 'rotated']
dpi = 3000

SaveData = False
PlotResults = False
Parallelization = (True, 3)

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
        'Ntraj': Ntraj,
        'r': r,
        'initial_J': initial_J,
        'EnergyPS' = EnergyPS,
        'type_traj': type_traj,
        'dpi': dpi,
        'ode_solver': 'DOP853',
        'Tol': [1e-10, 1e-10],
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'contour_levels': 50,
        'darkmode': darkmode})
###################################################################################################
