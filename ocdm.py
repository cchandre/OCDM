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
from scipy.optimize import root_scalar, fmin
from ocdm_modules import run_method
from ocdm_dict import dict_list, Parallelization, Omega
import multiprocessing

def run_case(dict):
	case = DiaMol(dict)
	run_method(case)

def main():
	print('\033[92m  Optical Centrifuge for Diatomic Molecules (here: Cl2)  \033[00m')
	if Parallelization[0]:
		if Parallelization[1] == 'all':
			num_cores = multiprocessing.cpu_count()
		else:
			num_cores = min(multiprocessing.cpu_count(), Parallelization[1])
		pool = multiprocessing.Pool(num_cores)
		pool.map(run_case, dict_list)
	else:
		for dict in dict_list:
			run_case(dict)

class DiaMol:
	def __repr__(self):
		return '{self.__class__.__name__}({self.DictParams})'.format(self=self)

	def __str__(self):
		return 'Optical Centrifuge for Diatomic Molecules ({}) with E0 = {:.3e}'.format(self.__class__.__name__, self.E0)

	def __init__(self, dict):
		for key in dict:
			setattr(self, key, dict[key])
		self.DictParams = dict
		self.Omega = lambda t: Omega(t)
		t = sp.symbols('t')
		self.Phi = sp.lambdify(t, sp.integrate(self.Omega(t), t))
		self.te_ps = xp.asarray(self.te) / 2.418884254e-5
		a_s = [42.13, 15.4, 5.4, -5.25]
		a_m = [-1599.0948228665286, 1064.701691434201, -262.7958617988855, 31.287242627165202, -1.8164825476900417, 0.04141328363593082]
		r_a = [5, 10]
		b_s = [25.29, 2.87, -0.09, -0.42]
		b_m = [68.28948313113857, -56.39145030911223, 25.523842390023567, -5.33917063891859, 0.5409284377558191, -0.02155141584655734]
		r_b = [3, 6]
		self.re, self.De, self.gam = 3.756, 0.0915, 1.0755
		self.we, self.Be = 2.5502e-3, 1.1098e-6
		al_Cl, al_Cl2 = 15.5421, 2 * 15.5421
		self.mu = 32548.53
		r = sp.Symbol('r')
		eps = self.De * (1 - sp.exp(-self.gam * (r - self.re)))**2 - self.De
		d_eps = sp.diff(eps, r)
		self.eps = sp.lambdify(r, eps)
		self.d_eps = sp.lambdify(r, d_eps)
		self.p_phi_ion = xp.floor(xp.sqrt(-self.mu * fmin(lambda r: -r**3 * self.d_eps(r), 5, full_output=True, disp=False, xtol=1e-10, ftol=1e-10)[1]))
		poly = lambda r, a: sum(c * r**k for k, c in enumerate(a))
		deriv = lambda fun: sp.lambdify(r, sp.diff(fun, r))
		para_s, perp_s = poly(r - self.re, a_s), poly(r - self.re, b_s)
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
		d_al_para = lambda r: xp.where(r<=r_a[0], d_para_s(r), xp.where(r<=r_a[1], d_para_m(r), d_para_l(r)))
		self.d_al_perp = lambda r: xp.where(r<=r_b[0], d_perp_s(r), xp.where(r<=r_b[1], d_perp_m(r), d_perp_l(r)))
		self.Dal = lambda r: self.al_para(r) - self.al_perp(r)
		self.d_Dal = lambda r: d_al_para(r) - self.d_al_perp(r)
		self.ZVS = lambda r, theta, phi, t: -self.mu * self.Omega(t)**2 * r**2 * xp.sin(theta)**2 / 2 + self.eps(r) - self.E0**2 * self.env(t)**2 / 4 * (self.Dal(r) * xp.sin(theta)**2 * xp.cos(phi)**2 + self.al_perp(r))

	def eqn_H(self, t, y_):
		Eeff = self.E0**2 * self.env(t)**2 / 4
		if self.dim == 2:
			r, phi, p_r, p_phi = xp.split(y_, 4)
			dr = p_r / self.mu
			dphi = p_phi / (self.mu * r**2) - self.Omega(t)
			dp_r = p_phi**2 / (self.mu * r**3) - self.d_eps(r) + Eeff * (self.d_Dal(r) * xp.cos(phi)**2 + self.d_al_perp(r))
			dp_phi = -Eeff * self.Dal(r) * xp.sin(2 * phi)
			return xp.concatenate((dr, dphi, dp_r, dp_phi), axis=None)
		elif self.dim == 3:
			r, theta, phi, p_r, p_theta, p_phi = xp.split(y_, 6)
			dr = p_r / self.mu
			dtheta = p_theta / (self.mu * r**2)
			dphi = p_phi / (self.mu * r**2 * xp.sin(theta)**2) - self.Omega(t)
			dp_r = p_theta**2 / (self.mu * r**3) + p_phi**2 / (self.mu * r**3 * xp.sin(theta)**2) - self.d_eps(r) + Eeff * (self.d_Dal(r) * xp.sin(theta)**2 * xp.cos(phi)**2 + self.d_al_perp(r))
			dp_theta = p_phi**2 * xp.cos(theta) / (self.mu * r**2 * xp.sin(theta)**3) + Eeff * self.Dal(r) * xp.sin(2 * theta) * xp.cos(phi)**2
			dp_phi = -Eeff * self.Dal(r) * xp.sin(theta)**2 * xp.sin(2 * phi)
			return xp.concatenate((dr, dtheta, dphi, dp_r, dp_theta, dp_phi), axis=None)

	def energy(self, t, y_, field=True):
		if field:
			Eeff = self.E0**2 * self.env(t)**2 / 4
		if self.dim == 2:
			if len(y_) > 2 * self.dim:
				r, phi, p_r, p_phi = xp.split(y_, 4)
			else:
				r, phi, p_r, p_phi = y_
			H = (p_r**2 + p_phi**2 / r**2) / (2 * self.mu) + self.eps(r)
			if field:
				H += - self.Omega(t) * p_phi - Eeff * (self.Dal(r) * xp.cos(phi)**2 + self.al_perp(r))
		elif self.dim == 3:
			if len(y_) > 2 * self.dim:
				r, theta, phi, p_r, p_theta, p_phi = xp.split(y_, 6)
			else:
				r, theta, phi, p_r, p_theta, p_phi = y_
			H = (p_r**2 + p_theta**2 / r**2 + p_phi**2 / (r**2 * xp.sin(theta)**2)) / (2 * self.mu) + self.eps(r)
			if field:
				H += - self.Omega(t) * p_phi - Eeff * (self.Dal(r) * xp.sin(theta)**2 * xp.cos(phi)**2 + self.al_perp(r))
		return H

	def env(self, t):
		te = xp.cumsum(self.te_ps)
		if self.envelope == 'sinus':
			return xp.where(t<=0, 0, xp.where(t<=te[0], xp.sin(xp.pi * t / (2 * te[0]))**2, xp.where(t<=te[1], 1, xp.where(t<=te[2], xp.sin(xp.pi * (te[2] - t) / (2 * self.te_ps[2]))**2, 0))))
		elif self.envelope == 'const':
			return 1
		elif self.envelope == 'trapez':
			return xp.where(t<=0, 0, xp.where(t<=te[0], t / te[0], xp.where(t<=te[1], 1, xp.where(t<=te[2], (te[2] - t) / self.te_ps[2], 0))))

	def initcond(self, N):
		if isinstance(self.initial_conditions[0], str):
			if self.initial_conditions[0] == 'microcanonical':
				r = self.generate_r(self.eps, self.initial_conditions[1], N, self.r)
				theta = xp.pi * xp.random.random((2, N))
				phi = 2 * xp.pi * xp.random.random((2, N)) - xp.pi
				P = xp.sqrt(2 * self.mu * (self.initial_conditions[1] - self.eps(r)))
				if self.dim == 2:
					p_r = P * xp.cos(phi[1])
					p_phi = P * xp.sin(phi[1]) * r
					return xp.concatenate((r, phi[0], p_r, p_phi), axis=None)
				elif self.dim == 3:
					p_r = P * xp.cos(phi[1]) * xp.sin(theta[1])
					p_theta = P * xp.sin(phi[1]) * xp.sin(theta[1]) * r
					p_phi = P * xp.cos(theta[1]) * r * xp.sin(theta[0])
					return xp.concatenate((r, theta[0], phi[0], p_r, p_theta, p_phi), axis=None)
			elif self.initial_conditions[0] == 'microcanonical_J':
				p0 = xp.sqrt(self.initial_conditions[2] * (self.initial_conditions[2] + 1))
				Energy0 = self.we * (self.initial_conditions[1] + 0.5) + self.Be * p0**2 - self.De
				phi = 2 * xp.pi * xp.random.random(N) - xp.pi
				p_phi = p0 * xp.ones(N)
				func = lambda r: p0**2 / (2 * self.mu * r**2) + self.eps(r)
				r = self.generate_r(func, Energy0, N, self.r)
				p_r = xp.sign(2 * xp.random.random(N) - 1) * xp.sqrt(2 * self.mu * (Energy0 - self.eps(r)) - p_phi**2 / r**2)
				if self.dim == 2:
					return xp.concatenate((r, phi, p_r, p_phi), axis=None)
				elif self.dim ==3:
					theta = xp.pi * xp.random.random(N)
					p_theta = xp.zeros(N)
					return xp.concatenate((r, theta, phi, p_r, p_theta, p_phi), axis=None)
		elif all([isinstance(item, float) or isinstance(item, int) for item in self.initial_conditions]):
			return xp.asarray(self.initial_conditions).flatten('F')
		else:
			raise ValueError('Wrong value for initial_conditions')

	def check_dissociation(self, y_):
		if self.dim == 2:
			dissociated = xp.zeros(len(y_)//4, dtype=bool)
			for _, (r, phi, p_r, p_phi) in enumerate(zip(*xp.split(y_, 4))):
				if xp.abs(p_phi) < self.p_phi_ion:
					Eval = self.energy(0, [r, phi, p_r, p_phi], field=False)
					if Eval > 0:
						rs = root_scalar(lambda r_: self.d_eps(r_) * r_**3 - p_phi**2 / self.mu, bracket=[4.9, 30], xtol=1e-10, rtol=1e-10, method='brentq').root
						Es = p_phi**2 / (2 * self.mu * rs**2) + self.eps(rs)
						if r > rs or Eval > Es:
							dissociated[_] = True
				else:
					dissociated[_] = True
			return dissociated
		elif self.dim == 3:
			r, theta, phi, p_r, p_theta, p_phi = xp.split(y_, 6)
			return (r > 20)

	def cart2sph(self, y_):
		if self.dim == 2:
			x, y, px, py = xp.split(y_, 4)
			r, phi = xp.hypot(x, y), xp.arctan2(y, x)
			p_r = (x * px + y * py) / r
			p_phi = x * py - y * px
			return xp.concatenate((r, phi, p_r, p_phi))
		elif self.dim == 3:
			x, y, z, px, py, pz = xp.split(y_, 6)
			xy, phi = xp.hypot(x, y), xp.arctan2(y, x)
			r, theta = xp.hypot(xy, z), xp.arctan2(xy, z)
			p_r = (x * px + y * py + z * pz) / r
			p_theta = ((x * px + y * py) * z - pz * xy**2) / xy
			p_phi = x * py - y * px
			return xp.concatenate((r, theta, phi, p_r, p_theta, p_phi))

	def sph2cart(self, y_):
		if self.dim == 2:
			r, phi, p_r, p_phi = xp.split(y_, 4)
			x, y = r * xp.cos(phi), r * xp.sin(phi)
			px = p_r * xp.cos(phi) - p_phi * xp.sin(phi) / r
			py = p_r * xp.sin(phi) + p_phi * xp.cos(phi) / r
			return xp.concatenate((x, y, px, py))
		elif self.dim == 3:
			r, theta, phi, p_r, p_theta, p_phi = xp.split(y_, 6)
			x, y, z = r * xp.sin(theta) * xp.cos(phi), r * xp.sin(theta) * xp.sin(phi), r * xp.cos(theta)
			px = p_r * xp.sin(theta) * xp.cos(phi) + p_theta * xp.cos(theta) * xp.cos(phi) / r - p_phi * xp.sin(phi) / (r * xp.sin(theta))
			py = p_r * xp.sin(theta) * xp.sin(phi) + p_theta * xp.cos(theta) * xp.sin(phi) / r + p_phi * xp.cos(phi) / (r * xp.sin(theta))
			pz = p_r * xp.cos(theta) - p_theta * xp.sin(theta) / r
			return xp.concatenate((x, y, z, px, py, pz))

	def rotating(self, y_, angle, type='cartesian'):
		y__ = y_.copy()
		if type == 'spherical':
			y__  = self.sph2cart(y__)
		if self.dim == 2:
			x, y, px, py = xp.split(y__, 4)
		elif self.dim == 3:
			x, y, z, px, py, pz = xp.split(y__, 6)
		xr = x * xp.cos(angle) - y * xp.sin(angle)
		yr = x * xp.sin(angle) + y * xp.cos(angle)
		pxr = px * xp.cos(angle) - py * xp.sin(angle)
		pyr = px * xp.sin(angle) + py * xp.cos(angle)
		if self.dim == 2:
			y__ = xp.concatenate((xr, yr, pxr, pyr))
		elif self.dim == 3:
			y__ = xp.concatenate((xr, yr, z, pxr, pyr, pz))
		if type == 'spherical':
			return self.cart2sph(y__)
		return y__

	def mod(self, y_):
		return self.cart2sph(self.sph2cart(y_))

	def generate_r(self, func, Energy0, N, r):
		if func(xp.linspace(r[0], r[1], 2**12)).min() > Energy0:
			raise ValueError('Empty energy surface')
		else:
			vec = []
			while len(vec) <= N:
				vec_t = (r[1] - r[0]) * xp.random.random(N) + r[0]
				vec = xp.hstack((vec, vec_t[func(vec_t)<=Energy0]))
			return vec[:N]

if __name__ == "__main__":
	main()
