import os as os
import random
import string
import sys
from os import chdir, environ, getcwd, mkdir, path, sep, walk
from shutil import copy, copytree

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import openturns as ot
import openturns.viewer as viewer
import scipy as scipy
import scipy.stats as stats
import xarray as xr
from build_PCE import *

#########Load the data############
# stochastic inputs
stochastic_varnames = ["wind_direction", "wind_speed", "z0"]
plot_labels = ["dir", "speed", "z0"]

marginals = ["normal", "normal", "uniform"]
MC_sample_size = 100
number_of_turbines = 4  # TODO:read

power_file = "output_data_examples" + sep + "turbine_data_Stochastic_atHubheight.nc"
input_nc_file = nc.Dataset(
    "../cases/windio_4turbines/plant_energy_resource/Stochastic_atHubHeight.nc"
)
nvar = len(stochastic_varnames)
MC_sample = np.zeros((len(input_nc_file.variables["time"]), nvar))
for j in range(nvar):
    varname = stochastic_varnames[j]
    MC_sample[:, j] = input_nc_file.variables[varname][:]

# stochastic power outputs
power_data = nc.Dataset(power_file)
power_table = power_data.variables["power"][:, :]

#########Compute Sobol indices############
input_variable_array = MC_sample[:MC_sample_size, :]
power_1st_sobol_indices = np.zeros((number_of_turbines, nvar))
power_total_sobol_indices = np.zeros((number_of_turbines, nvar))
sample_std = np.std(power_table, axis=0)
PCE_deg = 3
#
copula_type = "independent"  # choices: gaussian, independent
for i in range(number_of_turbines):
    if number_of_turbines > 1:
        std_to_test = sample_std[i]
    else:
        std_to_test = sample_std
    if std_to_test != 0:
        copula = copula_type

        polynomialChaosResult = construct_PCE_ot(
            input_variable_array[:, :],
            power_table[i, :],
            marginals,
            copula,
            PCE_deg,
            LARS=True,
        )

        chaosSI = ot.FunctionalChaosSobolIndices(polynomialChaosResult)
        for v in range(nvar):
            if number_of_turbines > 1:
                power_1st_sobol_indices[i, v] = chaosSI.getSobolIndex(v)
                power_total_sobol_indices[i, v] = chaosSI.getSobolTotalIndex(v)
            else:
                power_1st_sobol_indices[i, v] = chaosSI.getSobolIndex(v)
                power_total_sobol_indices[i, v] = chaosSI.getSobolTotalIndex(v)


#########Plot power###############
fig, ax = plt.subplots(1, 4, figsize=(13, 3))
xs = np.linspace(0, 13, num=100)
for j in range(number_of_turbines):
    # plot and compare to histogram

    # plot and compare empirical distributions
    power_KDE = stats.gaussian_kde(power_table[:, j] / 10**6, bw_method="silverman")
    #
    ax[j].set_title("Turbine " + str(j + 1), fontsize=15)
    ax[j].grid(True)
    ax[j].plot(xs, power_KDE(xs))
    #
    ax[j].set_xlabel("Power (MW)")
    ax[j].set_ylabel("PDF")
    ax[j].tick_params(axis="both", which="major", labelsize=15)
    ax[j].xaxis.label.set_size(15)
    ax[j].yaxis.label.set_size(15)


plt.tight_layout()
plt.show()

#########Plot Sobol indices############
colors = ["#002d74", "#e85113"]
fig, ax = plt.subplots(1, 4, figsize=(13, 3))
# Bar positions
bar_width = 0.24
items = plot_labels
bar_positions1 = np.arange(len(items)) + 0.5 * bar_width
bar_positions2 = np.arange(len(items)) + 1.5 * bar_width

for i in range(number_of_turbines):
    values1 = []
    values2 = []
    for j in range(len(items)):
        values1.append(power_1st_sobol_indices[i, j])
        values2.append(power_total_sobol_indices[i, j])

    # Create bar plots
    ax[i].bar(
        bar_positions1,
        values1,
        width=bar_width,
        label="$1^{st}$ indices",
        color=colors[0],
    )
    ax[i].bar(
        bar_positions2, values2, width=bar_width, label="Total indices", color=colors[1]
    )

    # Set labels and title
    ax[i].set_xlabel("Variables")
    if i == 0:
        ax[i].set_ylabel(r"Sobol indices")
    # Set x-axis ticks and labels
    ax[i].set_xticks(bar_positions2)
    ax[i].set_xticklabels(items)
    ax[i].set_title("Turbine " + str(i + 1), fontsize=15)
    if i == 3:
        ax[i].legend(loc="upper right", fontsize=12)
    ax[i].tick_params(axis="both", which="major", labelsize=15)
    ax[i].xaxis.label.set_size(15)
    ax[i].yaxis.label.set_size(15)
    ax[i].set_ylim(0.0, 1.1)

# Show the plot
plt.tight_layout()
plt.show()
