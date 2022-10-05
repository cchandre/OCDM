# OCDM
Optical Centrifuge for Diatomic Molecules (OCDM) - codes for the chlorine molecule Cl<sub>2</sub>

- [`ocdm_dict.py`](https://github.com/cchandre/OCDM/blob/main/ocdm_dict.py): to be edited to change the parameters of the OCDM computation (see below for a dictionary of parameters)

- [`ocdm.py`](https://github.com/cchandre/OCDM/blob/main/ocdm.py): contains the DiaMol class and main functions defining the OCDM dynamics (for Cl<sub>2</sub>)

- [`ocdm_modules.py`](https://github.com/cchandre/OCDM/blob/main/ocdm_modules.py): contains the methods to integrate the OCDM dynamics

- [`read_DiaMol_dissocation.m`](https://github.com/cchandre/OCDM/blob/main/read_DiaMol_dissociation.m): MATLAB script to produce the dissociation probability figure from the output files `DiaMol_dissociation.txt` 

- [`read_DiaMol_trajectories.m`](https://github.com/cchandre/OCDM/blob/main/read_DiaMol_trajectories.m): MATLAB script to plot the trajectories from the output files `DiaMol.mat`

Once [`ocdm_dict.py`](https://github.com/cchandre/OCDM/blob/main/ocdm_dict.py) has been edited with the relevant parameters, run the file as 
```sh
python3 ocdm.py
```
or 
```sh
nohup python3 -u ocdm.py &>ocdm.out < /dev/null &
```
The list of Python packages and their version are specified in [`requirements.txt`](https://github.com/cchandre/OCDM/blob/main/requirements.txt)
___
##  Parameter dictionary

- *Method*: string
   - 'plot_potentials': plot the potential &epsilon;(*r*) and polarizabilities &alpha;(*r*)
   - 'plot_ZVS': plot zero velocity functions
   - 'dissociation': computes the dissociation probability as a function of the amplitude *F*<sub>0</sub>
   - 'trajectories': computes and displays the trajectories according to *type_traj*
   - 'poincaré'; Poincaré section is *&phi;*=0 (mod 2 &pi;) with *&phi;'*<0 in the plane (*r*,*p*<sub>*r*</sub>) if *EventPS*='phi', and *p<sub>*r*</sub> =0 with *p<sub>*r*</sub>'<0 in the plane (*&phi;*,*p*<sub>*&phi;*</sub>) if *EventPS*='phi'
- *dimension*: 2 or 3; dimension of the computation
- *F0*: float or array of floats; amplitude(s) *F*<sub>0</sub> of the electric field, *E*(*t*) = *F*<sub>0 </sub>*f*(*t*) [<b>e<sub>*x*</sub></b> cos&Phi;(*t*) + <b>e<sub>*y*</sub></b> sin&Phi;(*t*)] cos&omega;*t*, considered in the computation (atomic units)
- *Omega*: lambda function; values of the frequency of rotation of the polarisation axis, &Omega;=&Phi;'(*t*), as a function of time (atomic units)
- *envelope*: string ('const', 'sinus', 'trapez'); envelope function *f*(*t*) of the laser field
- *te*: array of 3 or 4 floats; duration of ramp-up, plateau, ramp-down and (optional) after pulse (in picoseconds)
- *Ntraj*: integer; number of trajectories to be integrated
- *r*: array of two floats; minimum and maximum values of *r* for the display of potentials, and range of *r* (atomic units) for the selection of initial conditions
- *initial_conditions*: string or array of floats; 
   - ['microcanonical', *E*<sub>0</sub>] for a microcanonical distribution with energy *E*<sub>0</sub>
   - ['microcanonical_J', *n*, *J*] for a microcanonical distribution with initial energy *E*<sub>0</sub> = &omega;<sub>e</sub> (*n*+1/2) + *B*<sub>e</sub> *J*(*J*+1)-*D*<sub>e</sub>
   - array of shape (*Ntraj*, 2*dimension*) containing the initial conditions to be integrated
- *spread3D*: float; between 0 and 1; spread in angle theta for the initial conditions (only in the 3D case)
- *EnergyPS*: float; initial value of the energy (atomic units) used in *Method*='poincaré'
- *EventPS*: string; 'phi' or 'pr'; choice of Poincaré section; Poincaré section is *&phi;*=0 (mod 2 &pi;) with *&phi;'*<0 in the plane (*r*,*p*<sub>*r*</sub>) if *EventPS*='phi', and *p<sub>*r*</sub> =0 with *p<sub>*r*</sub>'<0 in the plane (*&phi;*,*p*<sub>*&phi;*</sub>) if *EventPS*='phi'
- *ode_solver*: string; 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA', 'Verlet', 'BM4'; method for the integration of trajectories
    - 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA': (non-symplectic, variable time step); see [ivp_solve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) for more details
    - 'Verlet': Strang splitting or symplectic leap-frog integrator (symplectic, fixed time step, order 2); see [Wikipedia](https://en.wikipedia.org/wiki/Strang_splitting) for more details
    - 'BM4': BM<sub>6</sub>4 (symplectic, fixed time step, order 4) from [Blanes, Moan, J. Comput. Appl. Math. 142, 313 (2002)](https://doi.org/10.1016/S0377-0427(01)00492-7)
    - NB: For Poincaré sections, ode_solver = 'RK45' by default
- *ode_tol*: array of two floats; absolute and relative tolerances [atol, rtol] for variable time-step integrators; see [ivp_solve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) for more details
- *ode_step*: float; time step (in picoseconds) for the symplectic integrators 'Verlet' and 'BM4'
- *frame*: string; 'fixed' or 'rotating'; specifies in which frame the numerical integration is performed
- *type_traj*: array of 3 strings; ['all' or 'dissociated' or 'bounded', 'cartesian' or 'spherical', 'fixed' or 'rotating'] for the type of trajectories to be plotted and/or saved
- *dpi*: integer; dpi value for the figures 
- *criterion*: string; 'exact', 'distance' or 'angular_momentum'; criterion for dissociation
####
- *SaveData*: boolean; if True, the results are saved in a `.mat` file (with the type specified in 'type_traj' for 'trajectories'); NB: the dissociation probabilities are saved in a `.txt` file regardless of the value of SaveData
- *PlotResults*: boolean; if True, the results (for 'plot_potentials', 'plot_ZVS' and 'trajectories') are plotted right after the computation (with the type specified in 'type_traj' for 'trajectories')
- *Parallelization*: tuple (boolean, int); True for parallelization, int is the number of cores to be used or int='all' to use all available cores
- *darkmode*: boolean; if True, plots are done in dark mode

---
For more information: <cristel.chandre@cnrs.fr>

<p align="center">
  <img src="https://github.com/cchandre/OCDM/blob/main/Figure_2.png" alt="Example" width="400"/>
  <img src="https://github.com/cchandre/OCDM/blob/main/Figure_1.png" alt="Example" width="400"/>
</p>
