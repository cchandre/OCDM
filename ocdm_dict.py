###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp


Method = 'display_ZVS'
r = [30, 60]

dimension = 2

E0 = 5e-3
Omega = 0
envelope = 'const'
te = [30, 100, 30]

SaveData = False
PlotResults = True

dpi = 300
darkmode = True

###################################################################################################
##                             DO NOT EDIT BELOW                                                 ##
###################################################################################################
dict = {
        'Method': Method,
        'r': xp.linspace(r[0], r[1], dpi),
        'dim': dimension,
        'E0': E0,
        'Omega': Omega,
        'envelope': envelope,
        'te': xp.asarray(te) / 2.42e-5,
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'dpi': dpi,
        'darkmode': darkmode}
###################################################################################################
