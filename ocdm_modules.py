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
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.io import savemat
import time
from datetime import date
import os
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def run_method(case):
    if case.darkmode:
        cs = ['k', 'w', 'c', 'm', 'r']
    else:
        cs = ['w', 'k', 'b', 'm', 'r']
    if case.PlotResults:
        plt.rc('figure', facecolor=cs[0], titlesize=30, figsize=[8,8])
        plt.rc('text', usetex=True, color=cs[1])
        plt.rc('font', family='serif', size=24)
        plt.rc('axes', facecolor=cs[0], edgecolor=cs[1], labelsize=30, labelcolor=cs[1], titlecolor=cs[1])
        plt.rc('xtick', color=cs[1], labelcolor=cs[1])
        plt.rc('ytick', color=cs[1], labelcolor=cs[1])
        plt.rc('image', cmap='bwr')
    filestr = type(case).__name__
    if case.Method == 'plot_potentials' and case.PlotResults:
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
    elif case.Method == 'plot_ZVS' and case.PlotResults:
        filestr += '_' + 'E0{:.2e}'.format(case.E0).replace('.', '')
        fig, ax = plt.subplots(1, 1)
        r = xp.linspace(case.r[0], case.r[1], case.dpi)
        phi = xp.linspace(0, 2*xp.pi, case.dpi)
        Phi, R = xp.meshgrid(phi, r)
        ZVS = case.ZVS(R, xp.pi/2, Phi, 0)
        plt.contourf(Phi, R, ZVS, case.contour_levels, cmap=plt.cm.hot)
        plt.colorbar()
        plt.contour(Phi, R, ZVS, case.contour_levels, linewidths=1, colors='k', linestyles='solid')
        ax.set_xlabel(r'$\phi$')
        ax.set_ylabel(r'$r$')
        ax.set_title(r'Zero-Velocity Surface')
        if case.SaveData:
            fig.savefig(filestr + '.png', dpi=case.dpi)
            print('\033[90m        Figure saved in {}.png \033[00m'.format(filestr))
        plt.show()
    elif case.Method in ['dissociation', 'trajectories']:
        y0 = case.initcond(case.Ntraj)
        t_eval = xp.linspace(0, case.te.sum(), case.dpi)
        if case.Method == 'dissociation':
            t_eval = [t_eval[0], t_eval[-1]]
        if xp.any(y0):
            start = time.time()
            sol = solve_ivp(case.eqn_H, (t_eval[0], t_eval[-1]), y0, method=case.ode_solver, t_eval=t_eval, atol=case.Tol[0], rtol=case.Tol[1])
            print('\033[90m        Computation finished in {} seconds \033[00m'.format(int(time.time() - start)))
            dissociated = case.check_dissociation(sol.y[:, -1])
            if case.type_traj[1] == 'cartesian':
                yc = case.sph2cart(sol.y)
            elif case.type_traj[1] == 'spherical':
                yc = case.mod(sol.y)
            if case.type_traj[2] == 'lab':
                yc = case.rotating(yc, case.Phi(t_eval), type=case.type_traj[1])
            if case.type_traj[0] == 'dissociated' and xp.any(dissociated):
                yc = yc[xp.tile(dissociated, 2 * case.dim), :]
            elif case.type_traj[0] == 'not_dissociated' and (not xp.all(dissociated)):
                yc = yc[xp.logical_not(xp.tile(dissociated, 2 * case.dim)), :]
            elif case.type_traj[0] != 'all' and (case.PlotResults or case.SaveData):
                print('\033[33m          Warning: All trajectories are being displayed and/or saved \033[00m')
            save_data(case, yc, filestr)
            if case.Method == 'dissociation':
                proba = dissociated.sum() / case.Ntraj
                print('\033[96m          for E0 = {:.3e}, dissociation probability = {:.3e} \033[00m'.format(case.E0, proba))
                vec_data = [case.E0, proba]
                file = open(type(case).__name__ + '_' + case.Method + '.txt', 'a')
                if os.path.getsize(file.name) == 0:
                    file.writelines('%   E0           proba \n')
                file.writelines(' '.join(['{:.6e}'.format(data) for data in vec_data]) + '\n')
                file.close()
            elif case.Method == 'trajectories' and case.PlotResults:
                fig = plt.figure(figsize=(8, 10))
                if case.type_traj[1] == 'cartesian':
                    axs = fig.add_gridspec(2, hspace=0.2).subplots(sharex=True)
                elif case.type_traj[1] == 'spherical':
                    axs = fig.add_gridspec(3, hspace=0.2).subplots(sharex=True)
                te = xp.cumsum(case.te)
                for ax in axs:
                    for t in te:
                        ax.axvline(x=t, lw=1, color=cs[1])
                    ax.set_xlim((t_eval[0], t_eval[-1]))
                    ax.set_xlabel(r'$t$')
                panels = xp.asarray((1,) * case.dim + (0,) * case.dim)
                colors = xp.tile(cs[2:2+case.dim], 2)
                if case.type_traj[1] == 'cartesian':
                    if case.dim == 2:
                        labels = [r'$x$', r'$y$', r'$p_x$', r'$p_y$']
                        ylabels = [r'$p_x$, $p_y$', r'$x$, $y$']
                    elif case.dim == 3:
                        labels = [r'$x$', r'$y$', r'$z$', r'$p_x$', r'$p_y$', r'$p_z$']
                        ylabels = [r'$p_x$, $p_y$, $p_z$', r'$x$, $y$, $z$']
                elif case.type_traj[1] == 'spherical':
                    if case.dim == 2:
                        labels = [r'$r$', r'$\phi$', r'$p_r$', r'$p_\phi$']
                        ylabels = [r'$p_r$, $p_\phi$', r'$r$', r'$\phi$']
                        panels[1] = 2
                    elif case.dim == 3:
                        labels = [r'$r$', r'$\theta$', r'$\phi$', r'$p_r$', r'$p_\theta$', r'$p_\phi$']
                        ylabels = [r'$p_r$, $p_\theta$, $p_\phi$', r'$r$', r'$\theta$, $\phi$']
                        panels[1:3] = 2
                for k, coord in enumerate(xp.split(yc, 2 * case.dim, axis=0)):
                    if panels[k] == 2:
                        axs[panels[k]].set_ylim((-xp.pi, xp.pi))
                        axs[panels[k]].set_yticks([-xp.pi, 0, xp.pi])
                        axs[panels[k]].set_yticklabels([r'$-\pi$', r'0', r'$\pi$'])
                    axs[panels[k]].plot(t_eval, coord.transpose(), colors[k], label=labels[k])
                for ylabel, ax in zip(ylabels, axs):
                    handles, labels = ax.get_legend_handles_labels()
                    by_label = dict(zip(labels, handles))
                    ax.legend(by_label.values(), by_label.keys(), loc='upper right', labelcolor='linecolor')
                    ax.set_ylabel(ylabel)
                plt.show()
    elif case.Method == 'poincaré':
        filestr += '_' + 'E0{:.2e}'.format(case.E0).replace('.', '')
        t = sp.symbols('t')
        if sp.diff(case.Omega(t), t) != 0:
            print('\033[33m          Warning: The frequency Omega is not constant \033[00m')
        Omega = case.Omega(0)
        def event_ps(t, y_):
            return (y_[1] + xp.pi) % (2 * xp.pi) - xp.pi
        event_ps.direction = -1
        pr_max = lambda r: xp.sqrt(2 * case.mu * (case.mu * r**2 * Omega**2 / 2 + case.E0**2 / 4 * case.al_para(r) - case.eps(r) + case.Energy0))
        rand = xp.random.random((2, case.Ntraj))
        r = (case.r[1] - case.r[0]) * rand[0] + case.r[0]
        p_r = pr_max(r) * (2 * rand[1] - 1)
        p_phi = case.mu * r**2 * Omega + event_ps.direction * r * xp.sqrt(pr_max(r)**2 - p_r**2)
        y0_ = xp.vstack((r, xp.zeros(case.Ntraj), p_r, p_phi))
        y_events = []
        start = time.time()
        for y0 in y0_.transpose():
            sol = solve_ivp(case.eqn_H, (0, case.te.sum()), y0, events=event_ps, method=case.ode_solver, atol=case.Tol[0], rtol=case.Tol[1])
            if xp.any(y_events):
                y_events = xp.vstack((y_events, sol.y_events[0]))
            else:
                y_events = sol.y_events[0]
        print('\033[90m        Computation finished in {} seconds \033[00m'.format(int(time.time() - start)))
        save_data(case, y_events, filestr)
        if case.PlotResults:
            fig, ax = plt.subplots(1, 1)
            r = xp.linspace(case.r[0], case.r[1], case.dpi)
            ax.plot(r, pr_max(r), 'r')
            ax.plot(r, -pr_max(r), 'r')
            ax.plot(y_events[1:, 0], y_events[1:, 2], cs[1] + '.')
            ax.set_xlim(case.r)
            ax.set_xlabel(r'$r$')
            ax.set_ylabel(r'$p_r$')
            plt.show()

def save_data(case, data, filestr, info=[]):
    if case.SaveData:
        mdic = case.DictParams.copy()
        mdic.update({'data': data, 'info': info})
        mdic.update({'date': date.today().strftime(" %B %d, %Y\n"), 'author': 'cristel.chandre@cnrs.fr'})
        savemat(filestr, mdic)
        print('\033[90m        Results saved in {}.mat \033[00m'.format(filestr))
