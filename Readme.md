# PALAGRAM Monograph

This repository contains the data used in the research paper:

> Gadal C., Schneider J., Bonamy C., Chauchat J., Dossmann Y., Kiesgen de Richter S., Mercier M. J., Naaim-Bouvet F., Rastello M., and Lacaze L. **Chapter 16: Particle-laden gravity currents: the lock-release slumping regime at the laboratory scale.** *Submitted to AGU Monograph.* 


## Repository organization

```
  palagram_monograph
  │
  └───data: data are stored here
  │   └───input_data: input data as sent by everyone
  │       └─── ... : NETCDF files
  │   └───output_data: processed data output by analysis.py (also contains input_data)
  │       └─── ... : NETCDF files
  └───analysis:
      └───analysis.py: analysis code, that reads input_data and writes output_data
  └───paper: contains source files for article
  │   └───figures: contains source figures
  │       └─── ... : PDF files
  │       └─── figure_scripts: contains figure scripts that reads data in data/output_data and writes figures in paper/figures
  │            └─── *.py : python scripts for figures
  │   └─── ... : various files (.tex, .bib, ...)
  │   └─── main.pdf : article preprint

  ```

## Data organization

The CSV file `dataset_summary.csv` offers a summary of all runs and corresponding experimental parameters, allowing for easier access to the data.

The folder `data/output_data` contains 287 netcdf4 files corresponding to each experimental run used in the paper. For each run, the structure of the NetCDF file is the following:
- attributes:
  - particle_type: particle type used (silica sand, glass beads, etc..)
  - label: filename
  - lab: lab where this run has been performed
  - run_oldID: Old filename, corresponding to the experimental notebook
  - author: author(s) that acquired this run
  - setup: setup used to acquire the data. See article.
  - dataset: Dataset classification of this run, See paper.

- dimensions(sizes): time(n)

- variables(dimensions):
  - At(): Atwood number
  - Fr(): Froude number (adi. initial current velocity)
  - H0(): initial heavy fluid height inside the lock
  - H_a(): ambient fluid height outside the lock
  - L0(): streamwise lock length
  - L_1(): streamwise tank length after the lock
  - Re(): Reynolds number
  - S(): Settling number
  - St(): Stokes number
  - T_a(): ambient temperature
  - T_f(): heavy fluid temperature inside the lock
  - W0(): crossstream lock width
  - a(): lock aspect ratio
  - alpha(): bottom slope
  - d(): particle diameter
  - gprime(): specific gravity
  - lamb(): adi. attenuation parameter
  - nu_a(): ambient viscosity
  - nu_f(): heavy fluid lock viscosity
  - phi(): initial particle volume fraction inside the lock
  - rho_a(): ambient fluid density
  - rho_c(): heavy fluid mix density inside the lock
  - rho_f(): 
  - rho_p(): particle density
  - t('time',): time vector
  - t0(): characteristic timescale, t0 = L0/u0
  - u0(): characteristic velocity scale, u0 = sqrt(gprime*H0)
  - vs(): particle Stokes velocity
  - x_front('time',): front position vector

Variables can sometimes possess the following attributes:
  - unit: corresponding unit
  - std: error(s) on the given quantity, calculated by error propagation from measurement uncertainties using the `uncertainties` module (https://pythonhosted.org/uncertainties/) in Python.
  - comments: comments on the given quantity (definition, formulas, etc ..)

## Related works

- Gadal, C., Mercier, M. J., Rastello, M., & Lacaze, L. (2023). Slumping regime in lock-release turbidity currents. *Journal of Fluid Mechanics*, 974, A4. [doi:10.1017/jfm.2023.762](https://doi.org/10.1017/jfm.2023.762) 

- Gadal, C., Mercier, M., Rastello, M., & Lacaze, L. (2023). Data used in 'Slumping regime in lock-release turbidity currents' [Data set]. In Journal of Fluid Mechanics (Vol. 974, p. A4). *Zenodo*. https://doi.org/10.5281/zenodo.10058946

- Schneider, J., Dossmann, Y., Farges, O. et al. Investigation of particle laden gravity currents using the light attenuation technique. *Exp Fluids*, 64, 23 (2023). [doi:10.1007/s00348-022-03562-y](https://doi.org/10.1007/s00348-022-03562-y)

- Chauchat, J., Cheng, Z., Nagel, T., Bonamy, C., and Hsu, T.-J. (2017) SedFoam-2.0: a 3-D two-phase flow numerical model for sediment transport, *Geosci. Model Dev.*, 10, 4367-4392, [doi:10.5194/gmd-10-4367-2017](https://doi.org/10.5194/gmd-10-4367-2017) and [github](https://github.com/sedfoam/sedfoam)
