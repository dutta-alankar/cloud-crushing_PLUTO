# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 14:59:59 2024

@author: alankar
"""
import numpy as np
from scipy.interpolate import interp1d
import os
import re
import subprocess as sp

from pygrackle import \
    chemistry_data, \
    evolve_constant_density, \
    setup_fluid_container

from pygrackle.utilities.model_tests import \
    get_model_set, \
    model_test_format_version

grackle_data_dir = "/freya/ptmp/mpa/adutt/pluto-mit-grackle/cloud-crushing_PLUTO/python-scripts/.venv/grackle/grackle_data_files/input"

## useful constants
yr = 365 * 24 * 60**2
Myr = 1e6 * yr
Gyr = 1e3 * Myr
pi = np.pi
pc = 3.0856775807e18
kpc = 1e3 * pc
Mpc = 1e3 * kpc
s = 1
cm = 1
K = 1
km = 1e5 * cm
mp = 1.67262192369e-24
me = 9.1093837e-28
mH = 1.6735e-24
kB = 1.3806505e-16
G = 6.6726e-8
H0 = 67.74
H0cgs = H0 * ((km / s) / Mpc)
dcrit0 = 3 * H0cgs**2 / (8.0 * pi * G)
MSun = 2.0e33

# Set the simulation parameters
eta  = 100
mach = 1.5
tcoolmBtcc = 0.50
Tcl = 1.0e+04 # K
cloud_pos = 10.0 # Rcl
metallicity = 1.0 # ZSun
gamma = 5/3.
ncl = 1.0e-03 # cm^-3

wind_extent = 100 # Rcl
prp_extent  = 22 # Rcl
RclBdcell = 16

redshift = 0.
specific_heating_rate = 0.
volumetric_heating_rate = 0.
# dictionary to store extra information in output dataset
extra_attrs = {}

# Set solver parameters
my_chemistry = chemistry_data()
my_chemistry.use_grackle = 1
my_chemistry.with_radiative_cooling = 0
my_chemistry.primordial_chemistry = 0
my_chemistry.metal_cooling = 1
my_chemistry.UVbackground = 1
my_chemistry.self_shielding_method = 0
my_chemistry.H2_self_shielding = 0
my_chemistry.grackle_data_file = os.path.join(grackle_data_dir, "CloudyData_UVB=HM2012.h5")

my_chemistry.use_specific_heating_rate = 1
my_chemistry.use_volumetric_heating_rate = 1

# Set units
my_chemistry.comoving_coordinates = 0 # proper units
my_chemistry.a_units = 1.0
my_chemistry.a_value = 1.0 / (1.0 + redshift) / \
my_chemistry.a_units
my_chemistry.density_units = mp
my_chemistry.length_units = kpc    
my_chemistry.time_units = Myr
my_chemistry.set_velocity_units()

# Call convenience function for setting up a fluid container.
# This container holds the solver parameters, units, and fields.
metal_mass_fraction = metallicity * my_chemistry.SolarMetalFractionByMass
temperature = np.logspace(1, 9, 200)
fc = setup_fluid_container(
    my_chemistry,
    density=mp,
    temperature=temperature,
    metal_mass_fraction=metal_mass_fraction,
    converge=True)

if my_chemistry.use_specific_heating_rate:
    fc["specific_heating_rate"][:] = specific_heating_rate
if my_chemistry.use_volumetric_heating_rate:
    fc["volumetric_heating_rate"][:] = volumetric_heating_rate

# get data arrays with symbolic units
data = fc.finalize_data()

LAMBDA = interp1d(data["temperature"], np.abs(data["cooling_rate"]), fill_value="extrapolate")
mu_val = interp1d(data["temperature"], data["mean_molecular_weight"], fill_value="extrapolate")
np.savetxt("grackle-gas_prop.txt", np.vstack(( data["temperature"].value, 
                                                np.abs(data["cooling_rate"].value), 
                                                data["mean_molecular_weight"].value )).T)

mu_wind = mu_val(eta*Tcl)
mu_cl = mu_val(Tcl)
chi = mu_cl/mu_wind * eta
print(f"mu_w = {mu_wind:.3f}, mu_cl = {mu_cl:.3f}")
print(f"chi = {chi:.2f}")

nwind = ncl/eta
Twind = eta*Tcl

vwind = mach*np.sqrt(gamma*kB*Twind/(mu_wind*mp))
Pwind = nwind*kB*Twind
Tmix  = np.sqrt(Tcl*Twind)
nmix  = np.sqrt(ncl*nwind)
Pmix  = nmix*kB*Tmix
nHmix = nmix*mu_val(Tmix)*(mp/mH)*0.716
tcoolmix = (1./(gamma-1))*Pmix/(nHmix*nHmix*LAMBDA(Tmix))
Rcl = vwind*tcoolmix/(np.sqrt(chi)*tcoolmBtcc)

UNIT_DENSITY = nwind*mu_wind*mp
UNIT_LENGTH = Rcl
UNIT_VELOCITY = vwind

print(f"UNIT_DENSITY  = {UNIT_DENSITY/(mu_wind*mp):.2e} cm^-3")
print(f"UNIT_LENGTH   = {UNIT_LENGTH/pc:.2e} pc")
print(f"UNIT_VELOCITY = {UNIT_VELOCITY/1.0e+05:.2e} km s^-1")

ini_content = f"""
[Grid]

X1-grid    1     0.00        {(wind_extent*RclBdcell):d}         u        {wind_extent:.2f}
X2-grid    1     {(-0.5*prp_extent):.2f}        {(prp_extent*RclBdcell):d}         u        {(0.5*prp_extent):.2f}
X3-grid    1     {(-0.5*prp_extent):.2f}        {(prp_extent*RclBdcell):d}         u        {(0.5*prp_extent):.2f}
"""

print(ini_content)

details = os.uname()
year = sp.getoutput('date +"%Y"')
user = sp.getoutput("echo $USER")
work_dir = sp.getoutput("cd .. && pwd")
sys_name = details[0]
node_name = details[1]
release = details[2]
arch = details[-1]
byte = sp.getoutput('lscpu | grep Endian').split()[-2].lower()
version = details[-2]
compiler_details = sp.getoutput('cat ../makefile | grep "ARCH" | grep "="').split()[-1]
mpi_compiler = sp.getoutput('cat ../PLUTO/Config/Linux.mpicc.defs | grep "CC"').split()[-1]
c_compiler = re.sub('[\W_]+', '', sp.getoutput(f"{mpi_compiler} --version").split()[0])
pluto_loc = "./PLUTO"
pluto_ver = sp.getoutput(f'cd .. && grep -r "PLUTO_VERSION" {pluto_loc}/Src/* | grep "define"').split()[-1][1:-1]

sysconf = f"""
USER           = {user}
WORKING_DIR    = {work_dir}
SYSTEM_NAME    = {sys_name}
NODE_NAME      = {node_name}
RELEASE        = {release}
ARCH           = {arch}
BYTE_ORDER     = {byte}
VERSION        = {version}
PLUTO_DIR      = {pluto_loc}
PLUTO_VERSION  = "{pluto_ver}"
C_COMPILER     = {c_compiler}
MPI_C_COMPILER = {mpi_compiler}
"""

with open("../sysconf.out", "w") as ascii:
    ascii.write(sysconf[1:])
