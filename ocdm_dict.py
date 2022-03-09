###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp


Method = 'display_V2D'
r = [3.6, 4]

dimension = 2

E0 = 2e-3
omega = 2e-5

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
        'omega': omega,
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'dpi': dpi,
        'darkmode': darkmode}
###################################################################################################
