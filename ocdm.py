#
# BSD 2-Clause License
#
# Copyright (c) 2021, Cristel Chandre
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
		a_m = [-18648.4, 13674.7, -3950.49, 564.588, -39.996, 1.12517]
		r_a = [5.5, 7.8]
		b_s = [25.29, 2.87, -0.09, -0.42]
		b_m = [79.9758, -64.2146, 26.4758, -4.99457, 0.448782, -0.0156147]
		r_b = [3.25, 7.25]
		re, De, gam, al_Cl, al_Cl2 = 3.757, 0.091, 0.566, 15.5421, 2 * 15.5421
		mu = 32545.85
		r = sp.Symbol('r')
		eps = De * (1 - sp.exp(-gam * (r - re)))**2 - De
		d_eps = sp.diff(eps, r)
		self.eps = sp.lambdify(r, eps)
		self.d_eps = sp.lambdify(r, d_eps)
		poly = lambda r, a: sum(c * r**k for k, c in enumerate(a))
		deriv = lambda a: sp.lambdify(r, sp.diff(a, r))
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
		self.V2D = lambda phi, r: -mu * self.omega**2 * r**2 / 2 + self.eps(r) - self.E0**2 / 4 * (self.Dal(r) * xp.cos(phi)**2 + self.al_perp(r))

	def eqn_H2D(self, t, y):
		r, phi, pr, pphi = xp.split(y, 4)
		dr = pr / self.mu
		dphi = pphi / (self.mu * r**2) - self.omega
		dpr = pphi**2 / (self.mu * r**3) - self.d_eps(r) + self.E0**2 / 4 * (self.d_Dal(r) * xp.cos(phi)**2 + self.d_al_perp(r))
		dpphi = -self.E0**2 / 4 * self.Dal(r) * xp.sin(2 * phi)
		return xp.concatenate((dr, dphi, dpr, dpphi), axis=None)

	def env(self, t):
		te = xp.cumsum(self.te)
		if self.envelope == 'sinus':
			return xp.where(t<=0, 0, xp.where(t<=te[0], xp.sin(xp.pi * t / (2 * te[0]))**2, xp.where(t<=te[1], 1, xp.where(t<=te[2], xp.sin(xp.pi * (t - te[2]) / (2 * self.te[2]))**2, 0))))
		elif self.envelope == 'const':
			return 1
		elif self.envelope == 'trapez':
			return xp.where(t<=0, 0, xp.where(t<=te[0], t / te[0], xp.where(t<=te[1], 1, xp.where(t<=te[2], (te[2] - t) / self.te[2], 0))))

if __name__ == "__main__":
	main()
