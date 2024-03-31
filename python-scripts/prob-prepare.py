# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 14:59:59 2024

@author: alankar
"""
import numpy as np
from scipy.interpolate import interp1d
import subprocess as sp
import os
import re

## useful constants
yr = 365 * 24 * 60**2
Myr = 1e6 * yr
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
X_solar = 0.7154
Y_solar = 0.2703
Z_solar = 0.0143

# Set the simulation features
cooling = True
catalyst = False
auto_compile = True
boost = True

# Set the simulation parameters
chi  = 100
eta  = 100
mach = 1.2
tcoolmBtcc = 0.5
Tcl = 1.0e+04 # K
cloud_pos = 4.0 # Rcl
metallicity = 1.0 # ZSun
gamma = 5/3.
ncl = 0.1 # cm^-3

wind_extent = 50 # Rcl
prp_extent  = 15 # Rcl
RclBdcell = 8

tcc = np.sqrt(chi)
dump_time = 1.0*tcc
analysis_time = 0.1*tcc
t_stop = 25*tcc

log_steps = 100

Xp = X_solar * (1 - metallicity * Z_solar) / (X_solar + Y_solar)
Yp = Y_solar * (1 - metallicity * Z_solar) / (X_solar + Y_solar)
Zp = metallicity * Z_solar

cooltable = np.loadtxt("../cooltable.dat")
LAMBDA = interp1d(cooltable[:,0], cooltable[:,1], fill_value="extrapolate")
mu = 1./(2*Xp+0.75*Yp+0.5625*Zp)
nwind = ncl/chi
Twind = eta*Tcl
vwind = mach*np.sqrt(gamma*kB*Twind/(mu*mp))
Pwind = nwind*kB*Twind
Tmix  = np.sqrt(Tcl*Twind)
nmix  = np.sqrt(ncl*nwind)
Pmix  = nmix*kB*Tmix
nHmix = nmix*mu*(mp/mH)*Xp
tcoolmix = (1./(gamma-1))*Pmix/(nHmix*nHmix*LAMBDA(Tmix))
Rcl = vwind*tcoolmix/(np.sqrt(chi)*tcoolmBtcc)

UNIT_DENSITY = nwind*mu*mp
UNIT_LENGTH = Rcl
UNIT_VELOCITY = vwind

print(f"UNIT_DENSITY  = {UNIT_DENSITY/(mu*mp):.2e} cm^-3")
print(f"UNIT_LENGTH   = {UNIT_LENGTH/pc:.2e} pc")
print(f"UNIT_VELOCITY = {UNIT_VELOCITY/1.0e+05:.2e} km s^-1")

output_dir = f"output-{'wb' if boost else 'nb'}-chi{chi:.1f}eta{eta:.1f}mach{mach:.2f}tcoolmBtcc{tcoolmBtcc:.2e}Tcl{Tcl:.2e}met{metallicity:.2f}-{'w_cool' if cooling else 'n_cool'}-res{RclBdcell}"

os.system(f"mkdir -p ../{output_dir}/Log_Files")
os.system(f"mkdir -p ../{output_dir}/snapshots")
if cooling:
    os.system("cp ../makefiles/makefile-tab_cool ../makefile")
else:
    os.system("cp ../makefiles/makefile-no_cool ../makefile")

def_content = f"""
#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        {'NO' if not(cooling) else 'TABULATED'}
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            7

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO
#define  INTERNAL_BOUNDARY              YES
#define  SHOW_TIMING                    NO
#define  SHOW_TIME_STEPS                YES
#define  BOOST                          {'YES' if boost else 'NO'} 

/* -- user-defined parameters (labels) -- */

#define  CHI                            0
#define  ETA                            1
#define  MACH                           2
#define  tCoolMbtCC                     3
#define  TCL                            4
#define  XOFFSET                        5
#define  ZMET                           6

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   {UNIT_DENSITY:.4e}
#define  UNIT_LENGTH                    {UNIT_LENGTH:.4e}
#define  UNIT_VELOCITY                  {UNIT_VELOCITY:.4e}

/* [End] user-defined constants (do not change this line) */
#define  MULTIPLE_LOG_FILES             YES
#define  VERBOSE                        NO
"""

with open("../definitions.h", "w") as ascii:
    ascii.write(def_content[1:])

ini_content = f"""
[Grid]

X1-grid    1     0.00        {(wind_extent*RclBdcell):d}         u        {wind_extent:.2f}
X2-grid    1     {(-0.5*prp_extent):.2f}        {(prp_extent*RclBdcell):d}         u        {(0.5*prp_extent):.2f}
X3-grid    1     {(-0.5*prp_extent):.2f}        {(prp_extent*RclBdcell):d}         u        {(0.5*prp_extent):.2f}

[Chombo Refinement]

Levels           4
Ref_ratio        2 2 2 2 2
Regrid_interval  2 2 2 2
Refine_thresh    0.3
Tag_buffer_size  3
Block_factor     8
Max_grid_size    64
Fill_ratio       0.75

[Time]

CFL              0.3
CFL_max_var      1.1
tstop            {t_stop:.1f}
first_dt         1e-07

[Solver]

Solver         hllc

[Boundary]

X1-beg        userdef
X1-end        outflow
X2-beg        outflow
X2-end        outflow
X3-beg        outflow
X3-end        outflow

[Static Grid Output]

uservar    5    Temp   ndens   PbykB   mach   cellvol 
output_dir ./snapshots
log_dir    ./Log_Files
dbl       -1.0          -1   single_file
flt       -1.0          -1   single_file
vtk       -1.0          -1   single_file
dbl.h5    {(10*dump_time):.2e}      -1   single_file
flt.h5    {dump_time:.2e}      -1   single_file
tab       -1.0          -1
ppm       -1.0          -1
png       -1.0          -1
log        {log_steps}
analysis  {analysis_time:.2e}      -1

[Chombo HDF5 output]

Checkpoint_interval  -1.0  0
Plot_interval         1.0  0

[Particles]

Nparticles             0   -1
particles_dbl       -1.0   -1
particles_flt       -1.0   -1
particles_vtk       -1.0   -1
particles_tab       -1.0   -1

[Parameters]

CHI                {chi:.1f}
ETA                {eta:.1f}
MACH               {mach:.2f}
tCoolMbtCC         {tcoolmBtcc:.2e}
TCL                {Tcl:.2e}
XOFFSET            {cloud_pos:.1f}
ZMET               {metallicity:.2f}
"""

with open("../pluto.ini", "w") as ascii:
    ascii.write(ini_content[1:])

details = sp.getoutput("uname -a").split()
year = sp.getoutput('date +"%Y"')
user = sp.getoutput("echo $USER")
work_dir = sp.getoutput("cd .. && pwd")
sys_name = details[0]
node_name = sp.getoutput("hostname")
release = details[2]
arch = details[-2]
byte = sp.getoutput('lscpu | grep Endian').split()[-2].lower()
version = f'{details[4]} {" ".join(details[6:11])} {year}'
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
PLUTO_VERSION  = {pluto_ver}
C_COMPILER     = {c_compiler}
MPI_C_COMPILER = {mpi_compiler}
"""

with open("../sysconf.out", "w") as ascii:
    ascii.write(sysconf[1:])

if catalyst:
    os.system("python generateCatalystAdaptor.py")

if auto_compile:
    os.system("cd .. && make -j8 && make clean")
    os.system(f"mv ../pluto.ini ../{output_dir}")
    os.system(f"mv ../definitions.h ../{output_dir}")
    os.system(f"mv ../sysconf.out ../{output_dir}")
    os.system(f"mv ../pluto ../{output_dir}")
    os.system(f"cp ../job-scripts/slurm-script ../{output_dir}")
    if cooling:
        os.system(f"cp ../cooltable.dat ../{output_dir}")

print(f"To run the job change directory using: \ncd ../{output_dir}")

