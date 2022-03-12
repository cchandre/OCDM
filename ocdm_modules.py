#
# BSD 2-Clause License
#
# Copyright (c) 2022, Cristel Chandre
# All rights reserved.
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
from scipy.integrate import solve_ivp
from scipy.io import savemat
import time
from datetime import date
import os

def run_method(case):
    if case.darkmode:
        cs = ['k', 'w', 'c', 'm', 'r']
    else:
        cs = ['w', 'k', 'b', 'm', 'r']
    plt.rc('figure', facecolor=cs[0], titlesize=30, figsize=[8,8])
    plt.rc('text', usetex=True, color=cs[1])
    plt.rc('font', family='serif', size=24)
    plt.rc('axes', facecolor=cs[0], edgecolor=cs[1], labelsize=30, labelcolor=cs[1], titlecolor=cs[1])
    plt.rc('xtick', color=cs[1], labelcolor=cs[1])
    plt.rc('ytick', color=cs[1], labelcolor=cs[1])
    plt.rc('image', cmap='bwr')
    filestr = type(case).__name__
    if case.Method == 'plot_potentials':
        plt.rcParams.update({'figure.figsize': [16, 8]})
        fig, axs = plt.subplots(1, 2)
        r = xp.linspace(case.r[0], case.r[1], case.dpi)
        axs[0].plot(r, case.eps(r), cs[1], lw=3, label=r'$\varepsilon(r)$')
        axs[1].plot(r, case.al_para(r), cs[2], lw=3, label=r'$\alpha_\parallel(r)$')
        axs[1].plot(r, case.al_perp(r), cs[3], lw=3, label=r'$\alpha_\perp(r)$')
        for ax in axs:
            ax.set_xlabel(r'$r$')
            ax.legend(loc='upper right', labelcolor='linecolor')
        if case.SaveData:
            fig.savefig(filestr + '.png', dpi=case.dpi)
            print('\033[90m        Figure saved in {}.png \033[00m'.format(filestr))
        plt.show()
    elif case.Method == 'plot_ZVS':
        filestr += '_' + 'E0{:.2e}'.format(case.E0).replace('.', '')
        fig, ax = plt.subplots(1, 1)
        r = xp.linspace(case.r[0], case.r[1], case.dpi)
        phi = xp.linspace(0, 2*xp.pi, case.dpi)
        Phi, R = xp.meshgrid(phi, r)
        ZVS = case.ZVS(R, xp.pi/2, Phi, 0)
        plt.contourf(Phi, R, ZVS, 50, cmap=plt.cm.hot)
        plt.colorbar()
        plt.contour(Phi, R, ZVS, 50, linewidths=1, colors='k', linestyles='solid')
        ax.set_xlabel(r'$\phi$')
        ax.set_ylabel(r'$r$')
        ax.set_title(r'Zero-Velocity Surface')
        if case.SaveData:
            fig.savefig(filestr + '.png', dpi=case.dpi)
            print('\033[90m        Figure saved in {}.png \033[00m'.format(filestr))
        plt.show()
    elif case.Method in ['dissociation', 'plot_trajectories']:
        y0 = case.initcond(case.Ntraj)
        t_eval = xp.linspace(0, case.te.sum(), case.dpi)
        if case.Method == 'dissociation':
            t_eval = [t_eval[0], t_eval[-1]]
		if y0:
        	sol = solve_ivp(case.eqn_H, (t_eval[0], t_eval[-1]), y0, t_eval=t_eval, atol=case.Tol, rtol=case.Tol)
			dissociated = case.check_dissociation(sol.y[:, -1])
        	if case.Method == 'dissociation':
            	proba = dissociated.sum() / case.Ntraj
            	print('\033[96m          for E0 = {:.3e}, dissociation probability = {:.3e} \033[00m'.format(case.E0, proba))
            	vec_data = [case.E0, proba]
            	file = open(type(case).__name__ + '_' + case.Method + '.txt', 'a')
            	if os.path.getsize(file.name) == 0:
                	file.writelines('%   E0           proba \n')
                	file.writelines(' '.join(['{:.6e}'.format(data) for data in vec_data]) + '\n')
                	file.close()
        	elif case.Method == 'plot_trajectories':
            	fig = plt.figure(figsize=(8, 10))
            	axs = fig.add_gridspec(2, hspace=0.2).subplots(sharex=True)
            	te = xp.cumsum(case.te)
            	for ax in axs:
                	for t in te:
                    	ax.axvline(x=t, lw=1, color=cs[1])
                	ax.set_xlim((t_eval[0], t_eval[-1]))
                	ax.set_xlabel(r'$t$')
            	yc = case.pol2cart(sol.y).transpose()
				if case.type_traj == 'dissociated' and xp.any(dissociated):
					yc = yc[:, dissociated]
				elif case.type_traj == 'not_dissociated' and (not xp.all(dissociated)):
					yc = yc[:, xp.logical_not(dissociated)]
				elif case.type_traj != 'all':
					print('\033[96m     Warning: All trajectories are displayed \033[00m')   
            	if case.dim == 2:
                	x, y, px, py = xp.split(yc, 4, axis=1)
                	axs[0].plot(t_eval, x, cs[2], label=r'$x$')
                	axs[0].plot(t_eval, y, cs[3], label=r'$y$')
                	axs[1].plot(t_eval, px, cs[2], label=r'$p_x$')
                	axs[1].plot(t_eval, py, cs[3], label=r'$p_y$')
                	axs[0].set_ylabel(r'$x$, $y$')
                	axs[1].set_ylabel(r'$p_x$, $p_y$')
            	elif case.dim == 3:
                	x, y, z, px, py, pz = xp.split(yc, 6, axis=1)
                	axs[0].plot(t_eval, x, cs[2], label=r'$x$')
                	axs[0].plot(t_eval, y, cs[3], label=r'$y$')
                	axs[0].plot(t_eval, z, cs[4], label=r'$z$')
                	axs[1].plot(t_eval, px, cs[2], label=r'$p_x$')
                	axs[1].plot(t_eval, py, cs[3], label=r'$p_y$')
                	axs[1].plot(t_eval, pz, cs[4], label=r'$p_z$')
                	axs[0].set_ylabel(r'$x$, $y$, $z$')
                	axs[1].set_ylabel(r'$p_x$, $p_y$, $p_z$')
            	for ax in axs:
                	ax.legend(loc='upper right', labelcolor='linecolor')
            	plt.show()

def save_data(case, data, filestr, info=[]):
    if case.SaveData:
        mdic = case.DictParams.copy()
        mdic.update({'data': data, 'info': info})
        mdic.update({'date': date.today().strftime(" %B %d, %Y\n"), 'author': 'cristel.chandre@cnrs.fr'})
        savemat(filestr, mdic)
        print('\033[90m        Results saved in {}.mat \033[00m'.format(filestr))
