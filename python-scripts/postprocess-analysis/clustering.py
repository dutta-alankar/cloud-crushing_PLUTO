# -*- coding: utf-8 -*-
"""
Created on Fri May  10 14:16:08 2024

@author: alankar
"""
import numpy as np
import seaborn as sns
import os
import hdbscan
import h5py
import matplotlib.pyplot as plt
from decimal import Decimal

sns.set_context('poster')
sns.set_style('white')
sns.set_color_codes()
plot_kwds = {'alpha' : 0.5, 's' : 80, 'linewidths':0}
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
# fig = plt.figure(figsize=(13, 10))
# colors = sns.color_palette("Paired", n_colors=len(tcoolmBtcc))  # a list of RGB tuples

till = 25

for i, tcBtcc in enumerate(tcoolmBtcc):
    output_dir = f"../../output-{'wb' if boost else 'nb'}-chi{chi:.1f}eta{eta:.1f}mach{mach:.2f}tcoolmBtcc{tcBtcc:.2e}Tcl{Tcl:.2e}met{metallicity:.2f}-{'w_cool' if cooling else 'n_cool'}-res{RclBdcell}"
    if output_dir in run_dirs:
        print(output_dir)
        clustering_analysis = f"{output_dir}/clustering"
        os.makedirs(clustering_analysis, exist_ok = True)
        for file_no in range(till+1):
            with h5py.File(f"{output_dir}/snapshots/data.{file_no:04d}.flt.h5", "r") as data:
                temperature = np.array(data[f"/Timestep_{file_no}/vars/Temp"]).flatten()
                x_cells = np.array(data["/cell_coords/X"]).flatten()
                y_cells = np.array(data["/cell_coords/Y"]).flatten()
                z_cells = np.array(data["/cell_coords/Z"]).flatten()

            condition = temperature<8e+04 # K
            temperature = temperature[condition]
            x_cells = x_cells[condition]
            y_cells = y_cells[condition]
            z_cells = z_cells[condition]
            coordinates = np.vstack( (x_cells, y_cells, z_cells) ).T
            print(file_no, coordinates.shape, end="\r")
            if coordinates.shape[0] == 0:
                print("No cloud exists!")
                continue

            clusterer = hdbscan.HDBSCAN(min_cluster_size=max(5, int(0.01*coordinates.shape[0])), 
                                        gen_min_span_tree=True)
            clusterer.fit(coordinates)

            clusterer.minimum_spanning_tree_.plot(edge_cmap='viridis', 
                                                  edge_alpha=0.6, 
                                                  node_size=80, 
                                                  edge_linewidth=2)
            plt.savefig(f"{clustering_analysis}/span-tree.{file_no:04d}.svg")
            plt.close()

            '''
            clusterer.single_linkage_tree_.plot(cmap='viridis', colorbar=True)
            plt.savefig(f"{clustering_analysis}/hierarchy-all.{file_no:04d}.svg")
            plt.close()
            '''

            clusterer.condensed_tree_.plot(select_clusters=True, selection_palette=sns.color_palette())
            plt.savefig(f"{clustering_analysis}/hierarchy-hdbscan.{file_no:04d}.svg")
            plt.close()

            cluster_label = clusterer.labels_
            cluster_prob  = clusterer.probabilities_
            np.save( f"{clustering_analysis}/cluster-labels.{file_no:04d}.npy", np.vstack([x_cells, y_cells, z_cells, cluster_label, cluster_prob, temperature]) )
