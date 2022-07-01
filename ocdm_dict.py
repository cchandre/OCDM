###################################################################################################
##               Dictionary of parameters: https://github.com/cchandre/OCDM                      ##
###################################################################################################

import numpy as xp

Method = 'dissociation'

dimension = 2

E0 = xp.linspace(0.01, 0.035, 200)
Omega = lambda t: 3e-10 * t
envelope = 'sinus'
te = [5, 40, 5, 5]

Ntraj = 50000
r = [2.5, 10]
initial_conditions = ['microcanonical_J', 0, 30]
EnergyPS = []
EventPS = 'phi'

type_traj = ['dissociated', 'cartesian', 'rotated']
dpi = 3000

SaveData = False
PlotResults = False
Parallelization = (True, 50)

darkmode = True

###################################################################################################
##                              DO NOT EDIT BELOW                                                ##
###################################################################################################
if Method == 'poincaré':
    dimension = 2
    envelope = 'const'
    te = [0, sum(te), 0]
    type_traj = ['all', 'spherical', 'rotated']
dict_list = [{'Method': Method} for _ in xp.atleast_1d(E0)]
if not isinstance(initial_conditions[0], str):
    initial_conditions = xp.asarray(initial_conditions)
    Ntraj = len(initial_conditions.flatten()) // (2*dimension)
for dict, E in zip(dict_list, xp.atleast_1d(E0)):
    dict.update({
        'dim': dimension,
        'E0': E,
        'envelope': envelope,
        'te': te,
        'Ntraj': Ntraj,
        'r': r,
        'initial_conditions': initial_conditions,
        'EnergyPS': EnergyPS,
        'event': EventPS,
        'type_traj': type_traj,
        'dpi': dpi,
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'darkmode': darkmode})
    dict.setdefault('contour_levels', 50)
    dict.setdefault('ode_solver', 'RK45')
    dict.setdefault('Tol', [1e-10, 1e-10])
###################################################################################################
