# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 09:28:14 2024

@author: alankar
"""
import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
from decimal import Decimal
from plot_style import *

def fexp(number):
    (sign, digits, exponent) = Decimal(number).as_tuple()
    return len(digits) + exponent - 1

def fman(number):
    return Decimal(number).scaleb(-fexp(number)).normalize()

run_dirs = [ f.path for f in os.scandir("../../") if f.is_dir() ]
run_dirs = [this_dir for this_dir in run_dirs if "output-" in this_dir]

# Set the simulation parameters
chi  = 100
eta  = 100
mach = 1.5
tcoolmBtcc = [0.01, 0.1, 0.5, 0.6, 1.0, 3.0, 5.0, 6.0, 7.0, 9.0, 10.0, 100.0] 
Tcl = 4.0e+04 # K
cloud_pos = 10.0 # Rcl
metallicity = 1.0 # ZSun
gamma = 5/3.
ncl = 0.1 # cm^-3

cooling = True
boost = True
RclBdcell = 8

tcc = np.sqrt(chi)

# create matplotlib figure
fig = plt.figure(figsize=(13, 10))
colors = sns.color_palette("Paired", n_colors=len(tcoolmBtcc))  # a list of RGB tuples

for i, tcBtcc in enumerate(tcoolmBtcc):
    output_dir = f"../../output-{'wb' if boost else 'nb'}-chi{chi:.1f}eta{eta:.1f}mach{mach:.2f}tcoolmBtcc{tcBtcc:.2e}Tcl{Tcl:.2e}met{metallicity:.2f}-{'w_cool' if cooling else 'n_cool'}-res{RclBdcell}"
    if output_dir in run_dirs:
        data = np.loadtxt(f"{output_dir}/snapshots/analysis.dat")
        label = r"$%d \times 10^{%d}$"%(np.round(float(fman(tcBtcc))), fexp(tcBtcc))
        '''
        if int(fexp(tcBtcc))==0:
            label = r"$%d$"%fman(tcBtcc)
        '''
        plt.plot( data[:,0]/tcc, data[:,5], color=colors[i],
                  label=label )

plt.xlabel(r"Time [$t_{\rm cc}$]")
plt.ylabel(r"Cold cloud mass $M(T<3.3\ T_{\rm cl})$ [$M_{\rm ini}$]")

plt.xlim(xmin=-0.1, xmax=25)
plt.ylim(ymin=2.0e-02, ymax=14.0) 

plt.legend(loc="best", prop={"size": 20}, ncol=4, title=r"$t_{\rm cool,mix}/t_{\rm cc}$",
           framealpha=0.5, shadow=False, fancybox=True)
plt.yscale("log")
plt.savefig("cloud-mass.png")
