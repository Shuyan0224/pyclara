# Start-to-End-S2E-simulation-for-CLARA-FEBE

## Overview

This repository contains a full **Start-to-End (S2E) simulation pipeline** for beam 
dynamics at the CLARA-FEBE beamline

The workflow connects conventional accelerator tracking codes with plasma simulation 
tools to study **beam-driven plasma wakefield acceleration (PWFA)** using realistic 
electron beam distributions.

---

### Tools used

- ASTRA — injector modelling (optional)
- Elegant — beam tracking through linac and beamline
- FBPIC — plasma wakefield simulation
- Python — data processing and analysis
- Simframe — ASTeC base model description

## Repository Structure

```text
clara-s2e/
│
├── Info/             # project overview and documentation
├── Injector/         # ASTRA injector simulations
├── Output/           # simulation results and data
├── PostInjector/     # Simframe and Elegant input files
├── src/              # python project source code 
└── README.md

---
```

## Installation

Clone the repository:

```bash
git clone -b elegant https://github.com/Shuyan0224/pyclara.git S2E
cd S2E
pip install -e .
```
**Note:** Elegant/Pelegant must be installed separately. See [Docker](#docker) for a pre-built environment.
 
## File Preparation
 
Run this once before your first simulation to split the lattice and create the input files:
 
```python
import pyclara
 
inputdir = '/path/to/elegant'
 
lte = pyclara.Converters.elegant_lte(filename=f'{inputdir}/FEBE.lte')
lte.load()
lte1, lte2 = lte.splitter('CLA-FEC1-SIM-FOCUS-01')
lte1.writer(f'{inputdir}/FEBE1.lte')
lte2.writer(f'{inputdir}/FEBE2.lte')
 
ele = pyclara.Converters.elegant_ele()
ele.load(f'{inputdir}/FEBE.ele')
ele.write(f'{inputdir}/FEBE1.ele')
ele.write(f'{inputdir}/FEBE2.ele')
```
 
## Running the S2E Pipeline
 
```python
import h5py
import pyclara
 
inputdir = '/path/to/elegant'
run_dir  = '/path/to/output'
 
elegant1 = pyclara.Simulation.TrackerElegant.Elegant_runner(inputdir, 'FEBE1', run_dir)
fbpic    = pyclara.Simulation.TrackerFBPIC.Fbpic_runner(run_dir)
elegant2 = pyclara.Simulation.TrackerElegant.Elegant_runner(inputdir, 'FEBE2', run_dir)
 
s = pyclara.Simulation.TrackerBeamline.TrackerBeamline()
s.set_input_particles('/path/to/elegant/FEBE_input.sdds')
s.add_tracker(elegant1)
s.add_tracker(fbpic)
s.add_tracker(elegant2)
 
final_particles = s.track(save_step=True)
```
 
## Elegant_runner Parameters
 
```python
e = pyclara.Simulation.TrackerElegant.Elegant_runner(inputdir, 'FEBE1', run_dir)
 
# input — auto-detects type, pass either an h5py.File or a .sdds path
e.set_input(h5py.File('input.h5', 'r'))   # from FBPIC output
e.set_input('path/to/input.sdds')          # from existing SDDS file
```
 
Override `e.cmd` to change the number of cores or Pelegant path:
 
**Docker (default):**
```python
# no changes needed — Docker paths are set by default
e.run()
```
 
**HPC or local:**
```python
e.cmd = [
    'mpirun',       # or full path — find with: which mpirun
    '-np', '32',    # number of cores
    'Pelegant',     # or full path — find with: which Pelegant
    f'{e.name}.ele'
]
e.run()
```
 
To find the paths on your system:
```bash
which Pelegant
which mpirun
```
 
## Fbpic_runner Input Parameters
 
Default parameters are set for the CLARA-FEBE plasma stage. Override before running:
 
```python
f = pyclara.Simulation.TrackerFBPIC.Fbpic_runner(run_dir)
 
# parameters
f.set_Sim_control(
    use_cuda = False,   # False = CPU , True = GPU 
    n_order  = -1,      # -1  = infinite order, most accurate, for single GPU or CPU
                        # 32 = finite order, faster, for multi-GPU MPI
    Nm       = 2        # 1 = cylindrically symmetric (fastest, ideal beam)
                        # 2 = captures dipole asymmetry (recommended for realistic beams)
                        # 3 = captures quadrupole asymmetry (strongly misaligned beams)
)

f.set_plasma_density(n0=1e22)   # n0 is background plasma density [m^-3]

f.get_lambda_p() #checking the plasma wavelength

f.set_Moving_Window(
    zmax = 1.5,   # window spans [-zmax, +zmax]*lambda_p longitudinally
    rmax = 1,     # radial extent in units of lambda_p
    Nz   = 512,   # number of grid cells in the longitudinal direction
    Nr   = 96     # number of grid cells in the radial direction
                   # more cells = finer resolution, slower simulation
)

f.set_sim_length(total_propagation = 2e-3)   
                    # total propagation distance through plasma [m]
                    # N_steps is computed automatically from this

f.set_plasma_size(
    p_zmin = 1,   # plasma start position in units of lambda_p
                   # beam is injected before this point and this need indisde moving window.
    p_nz   = 2,   # macro-particles per cell longitudinally
    p_nr   = 2,   # macro-particles per cell radially
    p_nt   = 4    # macro-particles per cell azimuthally
                   # more particles per cell = less noise, slower simulation
)

f.set_beam_charge(beam_charge = 250)   
                        # total beam charge [pC]
                        # used to compute macro-particle weights

 
# inject a Gaussian witness beam
f.set_input_Gaussian(
    sigma_z          = 1.4e-5,        # longitudinal RMS [m]
    sigma_r          = 1.85e-5,       # transverse RMS [m]
    n_emit           = 4.426719e-6,   # normalised emittance [m·rad]
    n_macroparticles = 262144,       # number of macro-particles
    injection_plane  = 1             # injection plane in units of lambda_p
)

 
# or inject a real beam from beam tracking code
f.set_input_h5(elegant_output.h5)
 # use pyclara.Converters to convert any lattice file to the required HDF5 format
 
# run with diagnostics
f.run(
    diag_period = 0,              # steps between diagnostic saves
                                   # 0 = save only the final step
                                   # e.g. 100 = save every 100 steps
    fieldtype   = ['E', 'rho'],    # options: 'E', 'B', 'rho', 'J'
    extra_species = None          # additional species e.g. {'witness': witness_beam}
)
 
output = f.get_output_h5()
```
 
## Converters
 
**Split and write lattice files:**
 
```python
lte = pyclara.Converters.elegant_lte(filename='FEBE.lte')
lte.load()
lte1, lte2 = lte.splitter('ELEMENT_NAME')
lte1.writer('FEBE1.lte')
```
 
**Load and modify .ele files:**
 
```python
ele = pyclara.Converters.elegant_ele()
ele.load('FEBE.ele')
ele.run_setup['lattice'] = 'FEBE1.lte'
ele.write('FEBE1.ele')
```
 
**Convert FBPIC output → Elegant input:**
 
```python
sdds_obj = pyclara.Converters.fbpic2sdds(h5py.File('output.h5', 'r'))
sdds_obj.save('beam_input.sdds')
```
 
**Convert Elegant output → particle dict:**
 
```python
d = pyclara.Converters.sdds2fbpic('FEBE_output.sdds')
# keys: x, y, z, xp, yp, p, px, py, dt
```
 
**Compute Twiss parameters from FBPIC output:**
 
```python
d = pyclara.Converters.fbpic2twiss(h5py.File('output.h5', 'r'))
# keys: p_central, s_start, beta_x, beta_y, alpha_x, alpha_y, eta_x, eta_y, etap_x, etap_y
```
 
## Coordinate Conventions
 
| | Elegant | FBPIC |
|---|---|---|
| Longitudinal | `t` [s] | `z` [m] |
| Divergence | `xp`, `yp` = px/p | `ux`, `uy` = px/(m_e·c) |
| Momentum | `p` = |βγ| | `uz` = pz/(m_e·c) |
 
All conversions are handled automatically by `fbpic2sdds` and `sdds2fbpic`.
 
## Docker
 
A pre-built image with Elegant, Pelegant, and FBPIC is available:
 
```bash
docker build -t shuyan/alma9-elegant --platform linux/amd64 -f alma9-elegant.txt .
 
docker run -v /path/to/your/files:/S2E -p 8888:8889 -ti shuyan/alma9-elegant:latest 
# Inside docker    
       pip install -e /S2E/path/to/pyclara
# if you want to use Jupyter lab 
       jupyter lab --allow-root --no-browser --ip=0.0.0.0 --port=8889 
           --IdentityProvider.token=''"
```

Open `http://127.0.0.1:8888/lab` in your browser.

## Reading elegant 

The output files can be read using python (you need to have
sdds installed, e.g. `pip install sdds`):

```bash
import sdds

d = sdds.load('FEBE.twi')

# making a standard optics plot 
subplot(2,1,1)
plot(d.getColumnValueList('s'),d.getColumnValueList('betax'))
plot(d.getColumnValueList('s'),d.getColumnValueList('betay'))
subplot(2,1,2)
plot(d.getColumnValueList('s'),d.getColumnValueList('etax'))
plot(d.getColumnValueList('s'),d.getColumnValueList('etay'))

# getting twiss parameters at a named element (e.g betax at CLA-FED-DIA-BPM-01-DRIFT-02
betax = array(d.getColumnValueList("betax"))[array(d.getColumnValueList('ElementName')) == 'CLA-FED-DIA-BPM-01-DRIFT-02']
print(betax)

```
The twi and sig files of the PostInjector optics ```FEBE.twi``` and ```FEBE.sig``` has been placed in the 
Output directory for convenience。

## Converting optics 

For example from elegant to ImpactX. Move to the PostInjector directory and run 

```python
import pyclara
import yaml

f = open("lattice.yaml")
d = yaml.safe_load(f)
l = pyclara.Converters.yaml2impactx(d)
```

## Authors

- Prof. Stewart Boogert — Director, Cockcroft Institute
- Shuyan Wen — University of Manchester / Cockcroft Institute
