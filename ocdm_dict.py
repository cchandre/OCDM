###################################################################################################
##            Dictionary of parameters: https://github.com/cchandre/OCDM                         ##
###################################################################################################

import numpy as xp


Method = 'display_potentials'
r_display = [3, 15]

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
		'r_display': xp.linspace(r_display[0], r_display[1], dpi),
        'E0': E0,
        'omega': omega,
        'SaveData': SaveData,
        'PlotResults': PlotResults,
        'dpi': dpi,
        'darkmode': darkmode}
###################################################################################################
