# OCDM
Optical Centrifuge for Diatomic Molecules (OCDM) - codes for the chlorine molecule Cl<sub>2</sub>

- [`ocdm_dict.py`](https://github.com/cchandre/OCDM/blob/main/ocdm_dict.py): to be edited to change the parameters of the OCDM computation (see below for a dictionary of parameters)

- [`ocdm.py`](https://github.com/cchandre/OCDM/blob/main/ocdm.py): contains the DiaMol class and main functions defining the OCDM dynamics (for Cl<sub>2</sub>)

- [`ocdm_modules.py`](https://github.com/cchandre/OCDM/blob/main/ocdm_modules.py): contains the methods to integrate the OCDM dynamics

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

- *Method*: string; 'plot_potentials', 'plot_ZVS', 'dissociation', 'trajectories'
####
- *dimension*: 2 or 3; dimension of the computation
- *E0*: float or array of floats; amplitude(s) *E*<sub>0</sub> of the electric field, *E*(*t*) = *E*<sub>0 </sub>*f*(*t*) [<b>e<sub>*x*</sub></b> cos&Phi;(*t*) + <b>e<sub>*y*</sub></b> sin&Phi;(*t*)] cos&omega;*t*, considered in the computation (atomic units)
- *Omega*: lambda function; values of the frequency of rotation of the polarisation axis, &Omega;=&Phi;'(*t*), as a function of time (atomic units)
- *envelope*: string ('const', 'sinus', 'trapez'); envelope function *f*(*t*) of the laser field
- *te*: array of 3 floats; duration of ramp-up, plateau and ramp-down (in picoseconds)
- *Energy0*: float (negative); value of the initial energy (atomic units)
- *r*: array of two floats; minimum and maximum values of *r* for the display of potentials, and range of *r* for the selection of initial conditions (atomic units)
- *Ntraj*: integer; number of trajectories to be integrated
- *plot_traj*: array of two strings; ['all' or 'dissociated' or 'non_dissociated', 'cartesian' or 'spherical'] for the type of trajectories to be plotted
- *dpi*: integer; dpi value for the figures 
####
- *SaveData*: boolean; if True, the results are saved in a `.mat` file; NB: the dissociation probabilities are saved in a `.txt` file regardless of the value of SaveData
- *PlotResults*: boolean; if True, the results (for 'plot_potentials', 'plot_ZVS' and 'trajectories') are plotted right after the computation
- *Parallelization*: tuple (boolean, int); True for parallelization, int is the number of cores to be used or int='all' to use all available cores
- *darkmode*: boolean; if True, plots are done in dark mode
####
These options may be changed in [`ocdm_dict.py`](https://github.com/cchandre/OCDM/blob/main/ocdm_dict.py) (see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html" target="_blank">'solve_ivp'</a>):
- *ode_solver*: integration method to use (default='RK45')
- *Tol*: absolute and relative tolerances (default=[1e-6, 1e-3])

---
For more information: <cristel.chandre@cnrs.fr>
