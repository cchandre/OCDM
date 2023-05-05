###################################################################################################
##               Dictionary of parameters: https://github.com/cchandre/OCDM                      ##
###################################################################################################

import numpy as xp

Method = 'poincaré'

dimension = 2

#F0 = xp.linspace(0.007, 0.03, 500)
F0 = 0.0121042
Omega = lambda t: 3e-10 * t
envelope = 'sinus'
te = [5, 300, 5]

Ntraj = 100
r = [3, 5]
initial_conditions = ['microcanonical_J', 0, 30]
spread3D = 0.1
TimePS = 47
EnergyPS = -0.177643
EventPS = 'phi'
ode_solver = 'BM4'
ode_tol = [1e-8, 1e-8]
ode_step = 1e-2
r_thresh = 20
frame = 'rotating'

type_traj = ['all', 'spherical', 'rotating']
dpi = 1024

SaveData = False
PlotResults = True
Parallelization = 4

darkmode = True

###################################################################################################
##                              DO NOT EDIT BELOW                                                ##
###################################################################################################
if Method == 'poincaré':
    dimension = 2
    envelope = 'const'
    OmegaPS = Omega(TimePS / 2.418884254e-5)
    Omega = lambda t: OmegaPS
    te = [0, sum(te), 0]
    type_traj = ['all', 'spherical', 'rotating']
    ode_solver = 'RK45'
dict_list = [{'Method': Method} for _ in xp.atleast_1d(F0)]
if not isinstance(initial_conditions[0], str):
    initial_conditions = xp.asarray(initial_conditions)
    Ntraj = len(initial_conditions.flatten()) // (2*dimension)
for dict, F in zip(dict_list, xp.atleast_1d(F0)):
    dict.update({
        'dim': dimension,
        'F0': F,
        'envelope': envelope,
        'te': te,
        'Ntraj': Ntraj,
        'r': r,
        'initial_conditions': initial_conditions,
        'spread3D': min(max(spread3D, 0), 1),
        'r_thresh': r_thresh,
        'EnergyPS': EnergyPS,
        'EventPS': EventPS,
        'type_traj': type_traj,
        'dpi': dpi if dpi >= 2 else 2,
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'darkmode': darkmode,
        'contour_levels': 50,
        'ode_solver': ode_solver,
        'Tol': ode_tol,
        'Step': ode_step,
        'frame': frame})
###################################################################################################
