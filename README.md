# OCDM
Optical Centrifuge for Diatomic Molecules

- [`ocdm_dict.py`](https://github.com/cchandre/OCDM/blob/main/ocdm_dict.py): to be edited to change the parameters of the OCDM computation (see below for a dictionary of parameters)

- [`ocdm.py`](https://github.com/cchandre/OCDM/blob/main/ocdm.py): contains the DiaMol class and main functions defining the OCDM dynamics

- [`ocdm_modules.py`](https://github.com/cchandre/OCDM/blob/main/ocdm_modules.py): contains the methods to integrate the OCDM dynamics

Once [`ocdm_dict.py`](https://github.com/cchandre/OCDM/blob/main/ocdm_dict.py) has been edited with the relevant parameters, run the file as 
```sh
python3 ocdm.py
```
or 
```sh
nohup python3 -u ocdm.py &>ocdm.out < /dev/null &
```

___
##  Parameter dictionary

- *Method*: string; 'display_potentials', 'display_ZVS', 'dissociation'
####
- *r*: array of two floats; minimum and maximum values of *r* for the display of the potentials and selection of initial conditions
- *dimension*: 2 or 3; dimension of the computation
- *E0*: float; amplitude of the electric field
- *Omega*: lambda function; values of the frequency of rotation of the polarisation axis as a function of time
- *envelope*: string ('const', 'sinus', 'trapez'); envelope function of the laser field
- *te*: array of 3 floats; duration of ramp-up, plateau and ramp-down (in picoseconds)
- *Energy0*: float (negative); value of the initial energy
- *Ntraj*: integer; number of trajectories to be integrated
- *Tol*: float; tolerance of the numerical integration 
####
- *SaveData*: boolean; if True, the results are saved in a `.mat` file
- *PlotResults*: boolean; if True, the results are plotted right after the computation
- *Parallelization*: tuple (boolean, int); True for parallelization, int is the number of cores to be used or int='all' to use all available cores
- *dpi*: integer; dpi value for the figures
- *darkmode*: boolean; if True, plots are done in dark mode

---
For more information: <cristel.chandre@cnrs.fr>
