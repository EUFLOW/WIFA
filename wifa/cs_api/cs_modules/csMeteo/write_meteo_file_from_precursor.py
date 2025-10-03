#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 16:22:29 2023

@author: j26483
"""

import argparse
import sys
from os import sep

import numpy as np


def temp2theta(temp, z_or_dz, P0, Pref=1000.0, g=9.81, Rair=287.0, Cp=1005.0):
    """
    computes the potential temperature as function of the thermodynamical
    temperature in a hydrostatic atmosphere
    input:
    temp: array of temperatures (K)
    z_or_dz: either the altitudes (m) of the profile
             or the differences between levels
    P0: pressure at the lower level (the z[0] level)
        in same units as Pref
    Pref: the 1000 hPa for the conventional definition of potential
          temperature
    g : gravity m/s**2
    Rair: dry air perfect gas constant J/kg/Kelvin
    Cp  : dry air heat capacity J/kg/Kelvin
    """
    RsCp = Rair / Cp
    theta0 = (Pref / P0) ** (RsCp) * temp
    gamma = g / Cp
    # print("g=%f Cp=%f gamma=%f"%(g,Cp,gamma,))
    psi = 1.0 / temp
    try:
        z = np.array(z_or_dz)
        dz = (
            z[..., 1:] - z[..., :-1]
        )  # Note perso : l'usage de l'ellipsis ici ne sert a rien
        psi[..., 1:] = np.cumsum(psi[..., :-1] * dz, axis=-1)
    except:  # Note perso : ce try except me parait suspect
        dz = np.array(z_or_dz)
        psi[..., 1:] = np.cumsum(psi[..., :-1] * dz, axis=-1)
    psi[..., 0] = 0.0
    # print(psi)
    theta = theta0 * np.exp(gamma * psi)
    return theta


if __name__ == "__main__":
    """
    Main function of script
    """

    parser = argparse.ArgumentParser(
        description="Script to generate a wind farm GMSH mesh using salome"
    )
    parser.add_argument("--resu_folder", dest="resu_folder", type=str)
    parser.add_argument("--meteo_filename", dest="meteo_filename", type=str)
    parser.add_argument("--hub_height", dest="hub_height", type=float)
    parser.add_argument("--hub_dir", dest="hub_dir", type=float)
    parser.add_argument("--Lmoinv", dest="Lmoinv", type=float)
    parser.add_argument("--precntmax", dest="precntmax", type=int)
    #
    args = parser.parse_args()
    RESU_folder = args.resu_folder
    meteo_file_name = args.meteo_filename
    hub_height = args.hub_height
    hub_dir = args.hub_dir
    Lmoinv = args.Lmoinv
    precntmax = args.precntmax
    #
    sea_level_pressure = 101325.0

    # get data
    data = np.genfromtxt(
        RESU_folder + sep + "profiles" + sep + "profile1.csv",
        delimiter=",",
        skip_header=1,
    )
    # get header
    with open(RESU_folder + sep + "profiles" + sep + "profile1.csv") as f:
        first_line = f.readline()
    var_name = first_line.replace(" ", "").replace("\n", "").split(",")

    T = data[:, var_name.index("RealTemp")]
    u = data[:, var_name.index("Velocity_x")]
    v = data[:, var_name.index("Velocity_y")]
    w = data[:, var_name.index("Velocity_z")]
    k = data[:, var_name.index("k")]
    eps = data[:, var_name.index("epsilon")]
    z = data[:, var_name.index("z")]

    pref = 100000
    theta = temp2theta(T + 273.15, z, sea_level_pressure, Pref=pref)

    # REORIENT PROFILE
    idx_zalpha = np.argmin(np.abs(z - hub_height))
    rotation_angle = np.radians(270.0 - hub_dir) - np.arctan2(
        v[idx_zalpha], u[idx_zalpha]
    )
    rotation_matrix = np.array(
        [
            [np.cos(rotation_angle), -np.sin(rotation_angle)],
            [np.sin(rotation_angle), np.cos(rotation_angle)],
        ]
    )
    u_rot, v_rot = np.dot(rotation_matrix, np.array([u, v]))

    # WRITE FILE
    # l'ordre etant: s, x, y, z, epsilon, k, RealTemp, Velocity_x, Velocity_y, Velocity_z
    meteo_file = open(meteo_file_name, "w")
    meteo_file.write("/ year, quantile, hour, minute, second of the profile:\n")
    meteo_file.write(" 2023,  1, 1,  0,  0\n")
    meteo_file.write("/ location of the meteo profile in the domaine (x,y):\n")
    meteo_file.write(" 0.0  0.0\n")
    meteo_file.write("/ Sea level pressure\n")
    meteo_file.write(" " + str(sea_level_pressure) + "\n")
    meteo_file.write(
        "/ Temperature profile: number of altitudes,(alt.,T in celcius,H in kg/kg ,Nc in N/cm**3)\n"
    )
    meteo_file.write(" " + str(len(z) + 1) + "\n")
    delimiter = "  "
    for i in range(len(z)):
        meteo_file.write(
            " "
            + str(z[i])
            + delimiter
            + str(T[i])
            + delimiter
            + str(0.0)
            + delimiter
            + str(0.0)
            + "\n"
        )
    # last line extrapolation to have h>1500m
    meteo_file.write(
        " "
        + str(z[-1] + 10)
        + delimiter
        + str(T[-1])
        + delimiter
        + str(0.0)
        + delimiter
        + str(0.0)
        + "\n"
    )

    meteo_file.write("/ Wind profile: number of altitudes,(alt.,u,v,ect,eps)\n")
    meteo_file.write(" " + str(len(z) + 1) + "\n")
    for i in range(len(z)):
        meteo_file.write(
            " "
            + str(z[i])
            + delimiter
            + str(u_rot[i])
            + delimiter
            + str(v_rot[i])
            + delimiter
            + str(k[i])
            + delimiter
            + str(eps[i])
            + "\n"
        )
    # last line extrapolation to have h>1500m
    meteo_file.write(
        " "
        + str(z[-1] + 10)
        + delimiter
        + str(u_rot[-1])
        + delimiter
        + str(v_rot[-1])
        + delimiter
        + str(k[-1])
        + delimiter
        + str(eps[-1])
        + "\n"
    )
    #
    meteo_file.close()
