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
        a_sr = [42.13, 15.4, 5.4, -5.25]
        a_mr = [-18648.4, 13674.7, -3950.49, 564.588, -39.996, 1.12517]
        r_a = [5.5, 7.8]
        b_sr = [25.29, 2.87, -0.09, -0.42]
        b_mr = [79.9758, -64.2146, 26.4758, -4.99457, 0.448782, -0.0156147]
        r_b = [3.25, 7.25]
        re, De, gam, al_Cl, al_Cl2 = 3.757, 0.091, 0.566, 15.5421, 2 * 15.5421
        r = sp.Symbol('r')
        eps = De * (1 - sp.exp(-gam * (r - re)))**2 - De
        d_eps = sp.diff(eps, r)
        self.eps = sp.lambdify(r, eps)
        self.d_eps = sp.lambdify(r, d_eps)
        poly = lambda r, a: sum(co*r**i for i,co in enumerate(a))
        deriv = lambda a: sp.lambdify(r, sp.diff(a, r))
        para_sr, perp_sr = poly(r - re, a_sr), poly(r - re, b_sr)
        para_mr, perp_mr = poly(r, a_mr), poly(r, b_mr)
        para_lr = (al_Cl2 + 4 * al_Cl**2 / r**3) / (1 - 4 * al_Cl**2 / r**6)
        perp_lr = (al_Cl2 - 2 * al_Cl**2 / r**3) / (1 - al_Cl**2 / r**6)
        d_para_sr, d_perp_sr = deriv(para_sr), deriv(perp_sr)
        d_para_mr, d_perp_mr = deriv(para_mr), deriv(perp_mr)
        d_para_lr, d_perp_lr = deriv(para_lr), deriv(perp_lr)
        para_sr, perp_sr = sp.lambdify(r, para_sr), sp.lambdify(r, perp_sr)
        para_mr, perp_mr = sp.lambdify(r, para_mr), sp.lambdify(r, perp_mr)
        para_lr, perp_lr = sp.lambdify(r, para_lr), sp.lambdify(r, perp_lr)
        self.al_para = lambda r: xp.where(r<=r_a[0], para_sr(r), xp.where(r<=r_a[1], para_mr(r), para_lr(r)))
        self.al_perp = lambda r: xp.where(r<=r_b[0], perp_sr(r), xp.where(r<=r_b[1], perp_mr(r), perp_lr(r)))
        self.d_al_para = lambda r: d_para_sr(r) if r <= r_a[0] else (d_para_mr(r) if r <= r_a[1] else d_para_lr(r))
        self.d_al_perp = lambda r: d_perp_sr(r) if r <= r_b[0] else (d_perp_mr(r) if r <= r_b[1] else d_perp_lr(r))
        self.Dal = lambda r: self.al_para(r) - self.al_perp(r)
        self.d_Dal = lambda r: self.d_al_para(r) - self.d_al_perp(r)

if __name__ == "__main__":
	main()
