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
import sympy as sp
from ocdm_modules import run_method
from ocdm_dict import dict

def run_case(dict):
	case = DiaMol(dict)
	run_method(case)

def main():
	run_case(dict)
	plt.show()

class DiaMol:
	def __repr__(self):
		return '{self.__class__.__name__}({self.DictParams})'.format(self=self)

	def __str__(self):
		return 'Optical Centrifuge for Diatomic Molecules ({self.__class__.__name__})'.format(self=self)

	def __init__(self, dict):
		for key in dict:
			setattr(self, key, dict[key])
		self.DictParams = dict
		a_s = [42.13, 15.4, 5.4, -5.25]
		a_m = [-18648.35842481222, 13674.730066523876, -3950.494622828384, 564.5884199329207, -39.995975029375124, 1.125168670489689]
		r_a = [5.5, 7.8]
		b_s = [25.29, 2.87, -0.09, -0.42]
		b_m = [79.97579598889766, -64.21463688187681, 26.47584281231729, -4.994574090432429, 0.4487817094919498, -0.015614680049778793]
		r_b = [3.75, 7.25]
		re, De, gam, al_Cl, al_Cl2 = 3.756775476853331, 0.0911, 0.566219624544, 15.5421, 2 * 15.5421
		self.mu = 32545.85
		r = sp.Symbol('r')
		eps = De * (1 - sp.exp(-gam * (r - re)))**2 - De
		d_eps = sp.diff(eps, r)
		self.eps = sp.lambdify(r, eps)
		self.d_eps = sp.lambdify(r, d_eps)
		poly = lambda r, a: sum(c * r**k for k, c in enumerate(a))
		deriv = lambda fun: sp.lambdify(r, sp.diff(fun, r))
		para_s, perp_s = poly(r - re, a_s), poly(r - re, b_s)
		para_m, perp_m = poly(r, a_m), poly(r, b_m)
		para_l = (al_Cl2 + 4 * al_Cl**2 / r**3) / (1 - 4 * al_Cl**2 / r**6)
		perp_l = (al_Cl2 - 2 * al_Cl**2 / r**3) / (1 - al_Cl**2 / r**6)
		d_para_s, d_perp_s = deriv(para_s), deriv(perp_s)
		d_para_m, d_perp_m = deriv(para_m), deriv(perp_m)
		d_para_l, d_perp_l = deriv(para_l), deriv(perp_l)
		para_s, perp_s = sp.lambdify(r, para_s), sp.lambdify(r, perp_s)
		para_m, perp_m = sp.lambdify(r, para_m), sp.lambdify(r, perp_m)
		para_l, perp_l = sp.lambdify(r, para_l), sp.lambdify(r, perp_l)
		self.al_para = lambda r: xp.where(r<=r_a[0], para_s(r), xp.where(r<=r_a[1], para_m(r), para_l(r)))
		self.al_perp = lambda r: xp.where(r<=r_b[0], perp_s(r), xp.where(r<=r_b[1], perp_m(r), perp_l(r)))
		self.d_al_para = lambda r: d_para_s(r) if r<=r_a[0] else (d_para_m(r) if r<=r_a[1] else d_para_l(r))
		self.d_al_perp = lambda r: d_perp_s(r) if r<=r_b[0] else (d_perp_m(r) if r<=r_b[1] else d_perp_l(r))
		self.Dal = lambda r: self.al_para(r) - self.al_perp(r)
		self.d_Dal = lambda r: self.d_al_para(r) - self.d_al_perp(r)
		self.ZVS = lambda r, theta, phi, t: -self.mu * self.Omega(t)**2 * r**2 * xp.sin(theta)**2 / 2 + self.eps(r) - self.E0**2 * self.env(t)**2 / 4 * (self.Dal(r) * xp.sin(theta)**2 * xp.cos(phi)**2 + self.al_perp(r))
		if -De < self.Energy0 < 0:
			self.rH0 = [re - xp.log(1 + xp.sqrt(1 + self.Energy0 / De)) / gam, re - xp.log(1 - xp.sqrt(1 + self.Energy0 / De)) / gam]
		else:
			self.rH0 = []

	def eqn_H2D(self, t, y):
		r, phi, pr, pphi = xp.split(y, 4)
		Eeff = E0**2 * self.env(t)**2 / 4
		dr = pr / self.mu
		dphi = pphi / (self.mu * r**2) - self.Omega(t)
		dpr = pphi**2 / (self.mu * r**3) - self.d_eps(r) + Eeff * (self.d_Dal(r) * xp.cos(phi)**2 + self.d_al_perp(r))
		dpphi = -Eeff * self.Dal(r) * xp.sin(2 * phi)
		return xp.concatenate((dr, dphi, dpr, dpphi), axis=None)

	def eqn_H3D(self, t, y):
		r, theta, phi, pr, ptheta, pphi = xp.split(y, 6)
		Eeff = E0**2 * self.env(t)**2 / 4
		dr = pr / self.mu
		dtheta = ptheta / (self.mu * r**2)
		dphi = pphi / (self.mu * r**2 * xp.sin(theta)**2) - self.Omega(t)
		dpr = ptheta**2 / (self.mu * r**3) + pphi**2 / (self.mu * r**3 * xp.sin(theta)**2) - self.d_eps(r) + Eeff * (self.d_Dal(r) * xp.cos(phi)**2 + self.d_al_perp(r))
		dptheta = pphi**2 * xp.cos(theta) / (self.mu * r**2 * xp.sin(theta)**3) + Eeff * self.Dal(r) * xp.sin(2 * theta) * xp.cos(phi)**2
		dpphi = -Eeff * self.Dal(r) * xp.sin(theta)**2 * xp.sin(2 * phi)
		return xp.concatenate((dr, dtheta, dphi, dpr, dptheta, dpphi), axis=None)

	def env(self, t):
		te = xp.cumsum(self.te)
		if self.envelope == 'sinus':
			return xp.where(t<=0, 0, xp.where(t<=te[0], xp.sin(xp.pi * t / (2 * te[0]))**2, xp.where(t<=te[1], 1, xp.where(t<=te[2], xp.sin(xp.pi * (te[2] - t) / (2 * self.te[2]))**2, 0))))
		elif self.envelope == 'const':
			return 1
		elif self.envelope == 'trapez':
			return xp.where(t<=0, 0, xp.where(t<=te[0], t / te[0], xp.where(t<=te[1], 1, xp.where(t<=te[2], (te[2] - t) / self.te[2], 0))))

	def initcond(self, N):
		if self.rH0:
			r0 = [max(self.r[0], self.rH0[0]), min(self.r[1], self.rH0[1])]
			r = (r0[1] - r0[0]) * xp.random.random(N) + r0[0]
			theta = xp.pi * xp.random.random(N)
			phi = 2 * xp.pi * xp.random.random(N)
			P = xp.sqrt(2 * self.mu * (self.Energy0 - self.eps(r)))
			Theta = xp.pi * xp.random.random(N)
			Phi = 2 * xp.pi * xp.random.random(N)
			if self.dimension == 2:
				pr = P * xp.cos(Phi)
				pphi = P * xp.sin(Phi)
				return xp.concatenate((r, phi, pr, pphi), axis=None)
			elif self.dimension == 3:
				pr = P * xp.cos(Phi) * xp.sin(Theta)
				ptheta = P * xp.sin(Phi) * xp.sin(Theta) * r
				pphi = P * xp.cos(Theta) * r * xp.sin(theta)
				return xp.concatenate((r, theta, phi, pr, ptheta, pphi), axis=None)
		else:
			print('\033[33m          Warning: Empty energy surface (H0>=0 or H0<=-De) \033[00m')
			return []



if __name__ == "__main__":
	main()
