#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as xp
import matplotlib.pyplot as plt

def run_method(case):
    if case.darkmode:
        cs = ['k', 'w', 'c', 'm']
    else:
        cs = ['w', 'k', 'b', 'm']
    plt.rc('figure', facecolor=cs[0], titlesize=30, figsize=[8,8])
    plt.rc('text', usetex=True, color=cs[1])
    plt.rc('font', family='serif', size=24)
    plt.rc('axes', facecolor=cs[0], edgecolor=cs[1], labelsize=30, labelcolor=cs[1], titlecolor=cs[1])
    plt.rc('xtick', color=cs[1], labelcolor=cs[1])
    plt.rc('ytick', color=cs[1], labelcolor=cs[1])
    plt.rc('image', cmap='bwr')
    print('\033[92m    {} \033[00m'.format(case.__str__()))
    print('\033[92m    E0 = {:.2f}   omega = {:.2f}  \033[00m'.format(case.E0, case.omega))
    filestr = type(case).__name__

    if case.Method == 'display_potentials':
        plt.rcParams.update({'figure.figsize': [16, 8]})
        fig, axs = plt.subplots(1, 2)
        axs[0].plot(case.r, case.eps(case.r), cs[1], lw=3, label=r'$\varepsilon(r)$')
        axs[1].plot(case.r, case.al_para(case.r), cs[2], lw=3, label=r'$\alpha_\parallel(r)$')
        axs[1].plot(case.r, case.al_perp(case.r), cs[3], lw=3, label=r'$\alpha_\perp(r)$')
        for ax in axs:
            ax.set_xlabel(r'$r$')
            ax.legend(loc='upper right', labelcolor='linecolor')
        if case.SaveData:
            fig.savefig(filestr + '.png', dpi=case.dpi)
            print('\033[90m        Figure saved in {}.png \033[00m'.format(filestr))
        plt.pause(0.5)
    elif case.Method == 'display_V2D':
        filestr += '_' + 'E0{:.2f}_OM{:.2f}'.format(case.E0, case.omega).replace('.', '')
        ig, ax = plt.subplots(1, 1)
        phi = xp.linspace(0, 2*xp.pi, 256)
        r = xp.linspace(case.r[0], case.r[1], 256)
        Phi, R = xp.meshgrid(phi, r)
        V2D = case.V2D(Phi, R)
        plt.contourf(Phi, R, V2D, label=r'$V_{2D}$')
        ax.set_xlabel(r'$\phi$')
        ax.set_ylabel(r'$r$')
        ax.legend(loc='upper right', labelcolor='linecolor')
        plt.pause(0.5)

#	elif case.Method == 'Poincare':

def save_data(case, data, filestr, info=[]):
	if case.SaveData:
		mdic = case.DictParams.copy()
		mdic.update({'data': data, 'info': info})
		mdic.update({'date': date.today().strftime(" %B %d, %Y\n"), 'author': 'cristel.chandre@cnrs.fr'})
		savemat(filestr, mdic)
		print('\033[90m        Results saved in {}.mat \033[00m'.format(filestr))
