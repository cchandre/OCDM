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
warnings.filterwarnings('ignore', category=RuntimeWarning)

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
    if case.Method == 'plot_potentials':
        plt.rcParams.update({'figure.figsize': [16, 8]})
        fig, axs = plt.subplots(1, 2)
        r = xp.linspace(case.r[0], case.r[1], case.dpi)
        axs[0].plot(r, case.eps(r), cs[1], lw=3, label=r'$\varepsilon(r)$')
        axs[1].plot(r, case.al_para(r), cs[2], lw=3, label=r'$\alpha_\parallel(r)$')
        axs[1].plot(r, case.al_perp(r), cs[3], lw=3, label=r'$\alpha_\perp(r)$')
        for ax in axs:
            ax.set_xlabel(r'$r$ (a.u.)')
            ax.legend(loc='upper right', labelcolor='linecolor')
            ax.set_xlim((case.r[0], case.r[1]))
        axs[0].set_ylim((-0.1, 0.2))
        if case.SaveData:
            fig.savefig(filestr + '.png', dpi=case.dpi)
            print(f'\033[90m        Figure saved in {filestr}.png \033[00m')
        plt.show()
    elif case.Method == 'plot_ZVS' and case.PlotResults:
        filestr += '_' + 'F0{:.2e}'.format(case.F0).replace('.', '')
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
            print(f'\033[90m        Figure saved in {filestr}.png \033[00m')
        plt.show()
    elif case.Method in ['dissociation', 'trajectories']:
        y0 = case.initcond(case.Ntraj)
        if xp.any(y0):
            start = time.time()
            if case.Method == 'dissociation':
                t_eval = xp.array([0, case.te_au.sum()])
            elif case.Method == 'trajectories':
                t_eval = xp.linspace(0, case.te_au.sum(), case.dpi)
            if case.ode_solver in ['RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA']:
                sol = solve_ivp(case.eqn_H, (t_eval[0], t_eval[-1]), y0, method=case.ode_solver, t_eval=t_eval, atol=case.Tol[0], rtol=case.Tol[1])
                yf = sol.y
            elif case.ode_solver in ['Verlet', 'BM4']:
                l = case.te_au.sum() / (case.dpi - 1)
                h = l / xp.ceil(l / case.Step)
                nstep = (case.dpi - 1) * int(xp.ceil(l / case.Step))
                eval = xp.zeros(nstep + 1, dtype=bool)
                if case.Method == 'dissociation':
                    eval[0], eval[-1] = True, True
                elif case.Method == 'trajectories':
                    eval[::int(xp.ceil(l / case.Step))] = True
                t, y = 0.0, y0.copy()
                yf = y0.copy()
                for _ in range(nstep):
                    t, y = case.eqn_sympl(h, t, y)
                    if eval[_ + 1]:
                        yf = xp.vstack((yf, y))
                yf = yf.transpose()
            print(f'\033[90m        Computation finished in {int(time.time() - start)} seconds \033[00m')
            dissociated, typeL = case.check_type(yf[:, -1])
            if case.type_traj[1] == 'cartesian':
                yc = case.sph2cart(yf)
            elif case.type_traj[1] == 'spherical':
                #yc = case.mod(yf)
                yc = yf.copy()
            if case.type_traj[2] == 'fixed' and case.frame == 'rotating':
                yc = case.rotating(yc, -case.Phi(t_eval), type=case.type_traj[1])
            elif case.type_traj[2] == 'rotating' and case.frame == 'fixed':
                yc = case.rotating(yc, case.Phi(t_eval), type=case.type_traj[1])
            if case.type_traj[0] == 'dissociated' and xp.any(dissociated):
                yc = yc[xp.tile(dissociated, 2 * case.dim), :]
            elif case.type_traj[0] == 'bounded' and (not xp.all(dissociated)):
                yc = yc[xp.logical_not(xp.tile(dissociated, 2 * case.dim)), :]
            elif case.type_traj[0] != 'all' and (case.PlotResults or case.SaveData):
                print('\033[33m          Warning: All trajectories are being displayed and/or saved \033[00m')
            t_eval *= 2.418884254e-5
            if case.type_traj[0] == 'all':
                save_data(case, xp.array([t_eval, yc[xp.tile(dissociated, 2 * case.dim), :], yc[xp.logical_not(xp.tile(dissociated, 2 * case.dim)), :]], dtype=object), filestr)
            else:
                save_data(case, xp.array([t_eval, yc], dtype=object), filestr)
            if case.Method == 'dissociation':
                diss_proba = dissociated.sum() / case.Ntraj
                typeL_proba = typeL.sum() / case.Ntraj
                print(f'\033[96m          for F0 = {case.F0:.3e}, dissociation probability = {diss_proba:.3e}, type-L probability = {typeL_proba:.3e} \033[00m')
                vec_data = [case.F0, diss_proba, typeL_proba]
                file = open(type(case).__name__ + '_' + case.Method + '.txt', 'a')
                if os.path.getsize(file.name) == 0:
                    file.writelines(f'%   initial = {case.initial_conditions}       beta = {case.beta:.3e}    dim = {case.dim}    N = {case.Ntraj}\n')
                    file.writelines(f'%   env = {case.envelope}     {case.te} \n')
                    file.writelines(f'%   F0           diss_proba     typeL_proba \n')
                file.writelines(' '.join([f'{data:.6e}' for data in vec_data]) + '\n')
                file.close()
            elif case.Method == 'trajectories' and case.PlotResults:
                fig = plt.figure(figsize=(12, 9.5))
                if case.type_traj[1] == 'cartesian':
                    axs = fig.add_gridspec(2, hspace=0.23).subplots(sharex=True)
                elif case.type_traj[1] == 'spherical':
                    axs = fig.add_gridspec(3, hspace=0.23).subplots(sharex=True)
                    axs[2].set_ylim((-1, 1))
                    axs[2].set_yticks([-1, 0, 1])
                    axs[2].set_yticklabels([r'-1', r'0', r'1'])
                te = xp.cumsum(xp.asarray(case.te))
                for ax in axs:
                    for t in te:
                        ax.axvline(x=t, lw=1, color=cs[1])
                    ax.set_xlim((t_eval[0], t_eval[-1]))
                    ax.set_xlabel(r'$t$ (ps)')
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
                        labels = [r'$r$', [r'cos $\phi$', r'sin $\phi$'], r'$p_r$', r'$p_\phi$']
                        ylabels = [r'$p_r$, $p_\phi$', r'$r$', r'cos $\phi$, sin $\phi$']
                        panels[1] = 2
                    elif case.dim == 3:
                        labels = [r'$r$', [r'cos $\theta$', r'sin $\theta$'], [r'cos $\phi$', r'sin $\phi$'], r'$p_r$', r'$p_\theta$', r'$p_\phi$']
                        ylabels = [r'$p_r$, $p_\theta$, $p_\phi$', r'$r$', r'cos $\theta$, sin $\theta$\\ cos $\phi$, sin $\phi$']
                        panels[1:3] = 2
                for k, coord in enumerate(xp.split(yc, 2 * case.dim, axis=0)):
                    if panels[k] == 2:
                        axs[2].plot(t_eval, xp.cos(coord.transpose()), colors[k], label=labels[k][0])
                        axs[2].plot(t_eval, xp.sin(coord.transpose()), colors[k], label=labels[k][1], linestyle='dashed')
                    else:
                        axs[panels[k]].plot(t_eval, coord.transpose(), colors[k], label=labels[k])
                r = xp.linspace(case.re, 4.25, case.dpi)
                axs[1].plot(xp.sqrt(case.mu * r**3 * case.d_eps(r)) / (case.mu * r**2 * case.beta) * 2.418884254e-5, r, cs[1])
                for ylabel, ax in zip(ylabels, axs):
                    handles, labels = ax.get_legend_handles_labels()
                    by_label = dict(zip(labels, handles))
                    ax.legend(by_label.values(), by_label.keys(), loc='upper right', labelcolor='linecolor')
                    ax.set_ylabel(ylabel, multialignment='center')
                plt.show()
    elif case.Method == 'poincaré':
        filestr += '_' + f'F0{case.F0:.2e}'.replace('.', '')
        t = sp.symbols('t')
        if sp.diff(case.Omega(t), t) != 0:
            print('\033[33m          Warning: The frequency Omega is not constant \033[00m')
        Omega = case.Omega(0)
        def event_ps(t, y_):
            if case.EventPS == 'phi':
                return (y_[1] + xp.pi) % (2 * xp.pi) - xp.pi
            elif case.EventPS == 'pr':
                return y_[2]
        event_ps.direction = -1
        if case.EventPS == 'phi':
            r = xp.linspace(case.r[0], case.r[1], 2**10)
            if case.EnergyPS < case.ZVS(r, xp.pi/2, 0, 0).min():
                raise ValueError('Empty Poincaré section')
            else:
                pr_max = lambda r: xp.sqrt(2 * case.mu * (case.EnergyPS - case.ZVS(r, 0, xp.pi/2, 0)))
                rand = xp.random.random((2, case.Ntraj))
                r = (case.r[1] - case.r[0]) * rand[0] + case.r[0]
                p_r = pr_max(r) * (2 * rand[1] - 1)
                p_phi = case.mu * r**2 * Omega + event_ps.direction * r * xp.sqrt(pr_max(r)**2 - p_r**2)
                y0_ = xp.vstack((r, xp.zeros(case.Ntraj), p_r, p_phi))
        elif case.EventPS == 'pr':
            r = xp.linspace(case.r[0], case.r[1], 2**10)
            phi = xp.linspace(0, 2 * xp.pi, 2**10)
            rp = xp.meshgrid(r, phi)
            if case.EnergyPS < case.ZVS(rp[0], xp.pi/2, rp[1], 0).min():
                raise ValueError('Empty Poincaré section')
            else:
                r, phi, p_phi = [], [], []
                while len(r) <= case.Ntraj:
                    rand = xp.random.random((2, case.Ntraj))
                    r_ = (case.r[1] - case.r[0]) * rand[0] + case.r[0]
                    phi_ = 2 * xp.pi * rand[1] - xp.pi
                    rp = xp.meshgrid(r_, phi_)
                    ZVS = case.ZVS(rp[0], xp.pi / 2, rp[1], 0)
                    r_, phi_, ZVS_ = r_[ZVS<=case.EnergyPS], phi_[ZVS<=case.EnergyPS], ZVS[ZVS<=case.EnergyPS]
                    r, phi = xp.hstack((r, r_)), xp.hstack((phi, phi_))
                    p_phi = case.mu * r_**2 * Omega + event_ps.direction * r_ * xp.sqrt(2 * case.mu * (case.EnergyPS - ZVS_))
                r, phi, p_phi = r[:case.Ntraj], phi[:case.Ntraj], p_phi[:case.Ntraj]
                y0_ = xp.vstack((r, phi, xp.zeros(case.Ntraj), p_phi))
        y_events = []
        start = time.time()
        for y0 in y0_.transpose():
            sol = solve_ivp(case.eqn_H, (0, case.te_au.sum()), y0, events=event_ps, method=case.ode_solver, atol=case.Tol[0], rtol=case.Tol[1])
            if xp.any(y_events):
                y_events = xp.vstack((y_events, sol.y_events[0]))
            else:
                y_events = sol.y_events[0]
        print(f'\033[90m        Computation finished in {int(time.time() - start)} seconds \033[00m')
        save_data(case, y_events, filestr)
        if case.PlotResults:
            fig, ax = plt.subplots(1, 1)
            if case.EventPS == 'phi':
                r = xp.linspace(case.r[0], case.r[1], case.dpi)
                ax.plot(r, pr_max(r), 'r')
                ax.plot(r, -pr_max(r), 'r')
                ax.plot(y_events[1:, 0], y_events[1:, 2], cs[1] + '.')
                ax.set_xlim(case.r)
                ax.set_xlabel(r'$r$')
                ax.set_ylabel(r'$p_r$')
            elif case.EventPS == 'pr':
                ax.plot((y_events[1:, 1] + xp.pi) %(2 * xp.pi) - xp.pi, y_events[1:, 3], cs[1] + '.')
                ax.set_xlim([-xp.pi, xp.pi])
                ax.set_xticks([-xp.pi, -xp.pi / 2, 0, xp.pi / 2, xp.pi])
                ax.set_xlabel(r'$\phi$')
                ax.set_ylabel(r'$p_\phi$')
                ax.set_xticklabels([r'$-\pi$', r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$'])
            plt.show()

def save_data(case, data, filestr, info=[]):
    if case.SaveData:
        mdic = case.DictParams.copy()
        mdic.update({'data': data, 'info': info})
        mdic.update({'date': date.today().strftime(' %B %d, %Y\n'), 'author': 'cristel.chandre@cnrs.fr'})
        savemat(filestr + '.mat', mdic)
        print(f'\033[90m        Results saved in {filestr}.mat \033[00m')
