#!/usr/bin/env python

###
### This file is generated for SALOME v9.10.0
###


import os
import sys

import salome

salome.salome_init()
from typing import List

import salome_notebook

notebook = salome_notebook.NoteBook()
import argparse
import math
from math import *

import GEOM
import numpy as np
import SALOMEDS
from salome.geom import geomBuilder

###
### GEOM component
###


geompy = geomBuilder.New()

parser = argparse.ArgumentParser(
    description="Script to generate a wind farm GMSH mesh using salome"
)
parser.add_argument("--wind_origin", dest="wind_origin", type=float)
parser.add_argument("--turbine_control", dest="turbine_control", type=float)
parser.add_argument("--output_file", dest="output_file", type=str)
parser.add_argument("--disk_mesh_size", dest="disk_mesh_size", type=float)
parser.add_argument("--domain_size", dest="domain_size", type=float)
parser.add_argument("--domain_height", dest="domain_height", type=float)
parser.add_argument("--damping_layer", dest="damping_layer", type=float)
args = parser.parse_args()

#######################
# lu dans des fichiers
#######################
xy_turbines = np.genfromtxt("Farm/DATA/turbines_info.csv", delimiter=",").transpose()[
    :2, :
]
nombre_turbines = xy_turbines.shape[1]

# dimensions en x et y du domaine
Lx = args.domain_size
Ly = args.domain_size

# caractéristiques turbines:
# hauteur moyeu
hub_heights_array = np.genfromtxt("Farm/DATA/turbines_info.csv", delimiter=",")[:, 2]
# diametre du rotor
diameters_array = np.genfromtxt("Farm/DATA/turbines_info.csv", delimiter=",")[:, 3]

# hm = np.max(hub_heights_array) #not used anymore, replaced by:
# Domain height
Lz = args.domain_height
extrusion1_max_height = np.max(hub_heights_array + diameters_array)
extrusion2_height = 0.6 * Lz
#
# angle de rotation (en degres)
provenance_vent = args.wind_origin
angle = 270.0 - provenance_vent

# nom du maillage généré
medFile = args.output_file

# choix domaine rectangulaire ou circulaire
dom = 0  # (0 pour circulaire et 1 pour rectangulaire)
###########################
# paramètres non utilisateur
###########################
# taille maille bord domaine
tc = args.domain_size / 150.0

# taille min
tm = args.disk_mesh_size
traff2 = 1.8 * tm
traff3 = 2.5 * tm
traff4 = 3.5 * tm
traff5 = 5.0 * tm
traff6 = 8.5 * tm

# epaisseur couche externe maillé en hexa
lb = 5 * tc
# caractéristiques maillage vertical
tv1 = 5
if args.damping_layer > 0:
    # TODO : generalize. Supposes a domain of height 25km
    tv2 = 2000
else:
    tv2 = 100
#########################
# début du script
#########################
origin = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(200, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 200, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 200)
geompy.addToStudy(origin, "origin")
geompy.addToStudy(OX, "OX")
geompy.addToStudy(OY, "OY")
geompy.addToStudy(OZ, "OZ")
Vz = geompy.MakeVectorDXDYDZ(0.0, 0.0, 1.0)
geompy.addToStudy(Vz, "Vz")

angle = np.radians(angle)

if dom == 0:
    domaine = geompy.MakeCylinderRH(max(Lx / 2, Ly / 2), Lz)
if dom == 1:
    domaine = geompy.MakeBoxDXDYDZ(Lx, Ly, Lz)
    geompy.TranslateDXDYDZ(domaine, -Lx / 2, -Ly / 2, 0)
geompy.addToStudy(domaine, "domaine")

part = geompy.MakeCylinderRH(max(Lx / 2, Ly / 2) - lb, Lz)

Point = geompy.MakeVertex(0, 0, 0)
Plane_1 = geompy.MakePlane(Point, OZ, Lx)

zone_raff1_eol = []
zone_raff2_eol = []
zone_raff3_eol = []
zone_raff4_pa = []
zone_raff5_pa = []
zone_raff6_pa = []
zone_raff7_pa = []

# placement des différentes zones de raffinement et zones sonde
if args.turbine_control > 0:
    for i in range(nombre_turbines):
        diametre = diameters_array[i]
        zone_raff1_eolienne = geompy.MakeBoxDXDYDZ(
            ((1.1 * diametre) // tm) * tm, ((1.1 * diametre) // tm) * tm, diametre + 40
        )
        S = geompy.TranslateDXDYDZ(
            zone_raff1_eolienne,
            -(((1.1 * diametre) // tm) * tm) / 2,
            -(((1.1 * diametre) // tm) * tm) / 2,
            0,
        )
        S = geompy.TranslateDXDYDZ(S, xy_turbines[0, i], xy_turbines[1, i], 0)
        zone_raff1_eol.append(S)

        zone_raff2_eolienne = geompy.MakeCylinderRH(
            ((1.1 * diametre) // tm) * tm + 4 * tm + 2 * traff2 * 0.5, diametre + 40
        )
        # S = geompy.TranslateDXDYDZ(zone_raff2_eolienne, -(8*tm+2*traff2)/4 , -(((1.1*diametre)//tm)*tm+4*tm+2*traff2)/4,0)
        S = geompy.MakeRotation(zone_raff2_eolienne, Vz, angle)
        S = geompy.TranslateDXDYDZ(S, xy_turbines[0, i], xy_turbines[1, i], 0)
        zone_raff2_eol.append(S)

        zone_raff3_eolienne = geompy.MakeCylinderRH(
            ((1.1 * diametre) // tm) * tm + 12 * tm + 2 * traff2 + 2 * traff3 * 0.5,
            diametre + 40,
        )
        # S = geompy.TranslateDXDYDZ(zone_raff3_eolienne, -(12*tm+2*traff2+2*traff3)/4, -(((1.1*diametre)//tm)*tm+12*tm+2*traff2+2*traff3)/4, 0)
        S = geompy.MakeRotation(zone_raff3_eolienne, Vz, angle)
        S = geompy.TranslateDXDYDZ(S, xy_turbines[0, i], xy_turbines[1, i], 0)
        zone_raff3_eol.append(S)

        zone_raff4_parc = geompy.MakeCylinderRH(8 * diametre * 0.5, diametre + 40)
        # S = geompy.TranslateDXDYDZ(zone_raff4_parc, -8*diametre/4, -8*diametre/4,0)
        S = geompy.TranslateDXDYDZ(
            zone_raff4_parc, xy_turbines[0, i], xy_turbines[1, i], 0
        )
        zone_raff4_pa.append(S)

        zone_raff5_parc = geompy.MakeCylinderRH(10.5 * diametre * 0.5, diametre + 40)
        # S = geompy.TranslateDXDYDZ(zone_raff5_parc, -10.5*diametre/4, -10.5*diametre/4,0)
        S = geompy.TranslateDXDYDZ(
            zone_raff5_parc, xy_turbines[0, i], xy_turbines[1, i], 0
        )
        zone_raff5_pa.append(S)

        zone_raff6_parc = geompy.MakeCylinderRH(
            (10 + 1.5 * tc / tm) * diametre * 0.5, diametre + 40
        )
        # S = geompy.TranslateDXDYDZ(zone_raff6_parc, -(10+1.5*tc/tm)*diametre/4, -(10+1.5*tc/tm)*diametre/4,0)
        S = geompy.TranslateDXDYDZ(
            zone_raff6_parc, xy_turbines[0, i], xy_turbines[1, i], 0
        )
        zone_raff6_pa.append(S)

        zone_raff7_parc = geompy.MakeCylinderRH(
            (10 + 3 * tc / tm) * diametre * 0.5, diametre + 40
        )
        # S = geompy.TranslateDXDYDZ(zone_raff7_parc, -(10+3*tc/tm)*diametre/4, -(10+3*tc/tm)*diametre/4,0)
        S = geompy.TranslateDXDYDZ(
            zone_raff7_parc, xy_turbines[0, i], xy_turbines[1, i], 0
        )
        zone_raff7_pa.append(S)
else:
    for i in range(nombre_turbines):
        diametre = diameters_array[i]
        zone_raff1_eolienne = geompy.MakeBoxDXDYDZ(
            4 * tm, ((1.1 * diametre) // tm) * tm, diametre + 40
        )
        S = geompy.TranslateDXDYDZ(
            zone_raff1_eolienne, -4 * tm / 4, -(((1.1 * diametre) // tm) * tm) / 4, 0
        )
        S = geompy.MakeRotation(zone_raff1_eolienne, Vz, angle)
        S = geompy.TranslateDXDYDZ(S, xy_turbines[0, i], xy_turbines[1, i], 0)
        zone_raff1_eol.append(S)

        zone_raff2_eolienne = geompy.MakeCylinderRH(
            ((1.1 * diametre) // tm) * tm + 4 * tm + 2 * traff2 * 0.5, diametre + 40
        )
        # S = geompy.TranslateDXDYDZ(zone_raff2_eolienne, -(8*tm+2*traff2)/4 , -(((1.1*diametre)//tm)*tm+4*tm+2*traff2)/4,0)
        S = geompy.MakeRotation(zone_raff2_eolienne, Vz, angle)
        S = geompy.TranslateDXDYDZ(S, xy_turbines[0, i], xy_turbines[1, i], 0)
        zone_raff2_eol.append(S)

        zone_raff3_eolienne = geompy.MakeCylinderRH(
            ((1.1 * diametre) // tm) * tm + 12 * tm + 2 * traff2 + 2 * traff3 * 0.5,
            diametre + 40,
        )
        # S = geompy.TranslateDXDYDZ(zone_raff3_eolienne, -(12*tm+2*traff2+2*traff3)/4, -(((1.1*diametre)//tm)*tm+12*tm+2*traff2+2*traff3)/4, 0)
        S = geompy.MakeRotation(zone_raff3_eolienne, Vz, angle)
        S = geompy.TranslateDXDYDZ(S, xy_turbines[0, i], xy_turbines[1, i], 0)
        zone_raff3_eol.append(S)

        zone_raff4_parc = geompy.MakeCylinderRH(8 * diametre * 0.5, diametre + 40)
        # S = geompy.TranslateDXDYDZ(zone_raff4_parc, -8*diametre/4, -8*diametre/4,0)
        S = geompy.TranslateDXDYDZ(
            zone_raff4_parc, xy_turbines[0, i], xy_turbines[1, i], 0
        )
        zone_raff4_pa.append(S)

        zone_raff5_parc = geompy.MakeCylinderRH(10.5 * diametre * 0.5, diametre + 40)
        # S = geompy.TranslateDXDYDZ(zone_raff5_parc, -10.5*diametre/4, -10.5*diametre/4,0)
        S = geompy.TranslateDXDYDZ(
            zone_raff5_parc, xy_turbines[0, i], xy_turbines[1, i], 0
        )
        zone_raff5_pa.append(S)

        zone_raff6_parc = geompy.MakeCylinderRH(
            (10 + 1.5 * tc / tm) * diametre * 0.5, diametre + 40
        )
        # S = geompy.TranslateDXDYDZ(zone_raff6_parc, -(10+1.5*tc/tm)*diametre/4, -(10+1.5*tc/tm)*diametre/4,0)
        S = geompy.TranslateDXDYDZ(
            zone_raff6_parc, xy_turbines[0, i], xy_turbines[1, i], 0
        )
        zone_raff6_pa.append(S)

        zone_raff7_parc = geompy.MakeCylinderRH(
            (10 + 3 * tc / tm) * diametre * 0.5, diametre + 40
        )
        # S = geompy.TranslateDXDYDZ(zone_raff7_parc, -(10+3*tc/tm)*diametre/4, -(10+3*tc/tm)*diametre/4,0)
        S = geompy.TranslateDXDYDZ(
            zone_raff7_parc, xy_turbines[0, i], xy_turbines[1, i], 0
        )
        zone_raff7_pa.append(S)

raff1_list = []
for i in range(nombre_turbines):
    raff1_list.append(zone_raff1_eol[i])
Union_raff1 = geompy.MakeFuseList(raff1_list, True, True)

raff2_list = []
for i in range(nombre_turbines):
    raff2_list.append(zone_raff2_eol[i])
Union_raff2 = geompy.MakeFuseList(raff2_list, True, True)

raff3_list = []
for i in range(nombre_turbines):
    raff3_list.append(zone_raff3_eol[i])
Union_raff3 = geompy.MakeFuseList(raff3_list, True, True)

raff4_list = []
for i in range(nombre_turbines):
    raff4_list.append(zone_raff4_pa[i])
Union_raff4 = geompy.MakeFuseList(raff4_list, True, True)

Union_raff5 = geompy.MakeFuseList(zone_raff5_pa, True, True)
Union_raff6 = geompy.MakeFuseList(zone_raff6_pa, True, True)
Union_raff7 = geompy.MakeFuseList(zone_raff7_pa, True, True)

if dom == 0:
    temp = [(i * max(Lx, Ly) / 2.0, 0.0) for i in [-1.0, 1.0]]
    temp = geompy.MakeEdge(
        geompy.MakeVertex(temp[0][0], temp[0][1], 0.0),
        geompy.MakeVertex(temp[1][0], temp[1][1], 0.0),
    )
    bande = geompy.MakeCut(domaine, part)
    part = geompy.MakePartition(
        [bande], [temp], [], [], geompy.ShapeType["FACE"], 0, [], 0
    )

part_list = []
part_list.append(Union_raff1)
part_list.append(Union_raff2)
part_list.append(Union_raff3)
part_list.append(Union_raff4)
part_list.append(Union_raff5)
part_list.append(Union_raff6)
part_list.append(Union_raff7)
if dom == 0:
    part_list.append(part)

Partition_2 = geompy.MakePartition([domaine], part_list)
Partition_2 = geompy.MakeCommonList([Partition_2, Plane_1], True)
geompy.addToStudy(Partition_2, "Partition_2")

# Creation des groupes (faces) pour CL
# Récuperation des paquets dedges en vis à vis

edges_x_raff1: List = []
edges_y_raff1: List = []
lListEdges = geompy.Propagate(Union_raff1)
for lEdges in lListEdges:
    edges = geompy.SubShapeAll(lEdges, geompy.ShapeType["EDGE"])
    vertex1 = geompy.MakeVertexOnCurve(edges[0], 0.25)
    vertex2 = geompy.MakeVertexOnCurve(edges[0], 0.75)
    basicCoords1 = geompy.PointCoordinates(vertex1)
    basicCoords2 = geompy.PointCoordinates(vertex2)
    tgtEdges = geompy.GetInPlace(Partition_2, lEdges)
    edges = geompy.SubShapeAll(tgtEdges, geompy.ShapeType["EDGE"])
    alongx_raff1 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
    alongy_raff1 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
    if (abs(basicCoords2[0] - basicCoords1[0])) > 5:
        # Selon X
        edges_x_raff1 = edges_x_raff1 + edges
    if (abs(basicCoords2[1] - basicCoords1[1])) > 5:
        # Selon Y
        edges_y_raff1 = edges_y_raff1 + edges
geompy.UnionList(alongx_raff1, edges_x_raff1)
geompy.UnionList(alongy_raff1, edges_y_raff1)

# Récuperation des paquets dedges en vis à vis
edges_x_raff2: List = []
edges_y_raff2: List = []
lListEdges = geompy.Propagate(Union_raff2)
for lEdges in lListEdges:
    edges = geompy.SubShapeAll(lEdges, geompy.ShapeType["EDGE"])
    vertex1 = geompy.MakeVertexOnCurve(edges[0], 0.25)
    vertex2 = geompy.MakeVertexOnCurve(edges[0], 0.75)
    basicCoords1 = geompy.PointCoordinates(vertex1)
    basicCoords2 = geompy.PointCoordinates(vertex2)
    tgtEdges = geompy.GetInPlace(Partition_2, lEdges)
    edges = geompy.SubShapeAll(tgtEdges, geompy.ShapeType["EDGE"])
    alongx_raff2 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
    alongy_raff2 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])

    if (abs(basicCoords2[0] - basicCoords1[0])) > 5:
        # Selon X
        edges_x_raff2 = edges_x_raff2 + edges
    if (abs(basicCoords2[1] - basicCoords1[1])) > 5:
        # Selon Y
        edges_y_raff2 = edges_y_raff2 + edges
geompy.UnionList(alongx_raff2, edges_x_raff2)
geompy.UnionList(alongy_raff2, edges_y_raff2)

# Récuperation des paquets dedges en vis à vis
edges_x_raff3: List = []
edges_y_raff3: List = []
lListEdges = geompy.Propagate(Union_raff3)
for lEdges in lListEdges:
    edges = geompy.SubShapeAll(lEdges, geompy.ShapeType["EDGE"])
    vertex1 = geompy.MakeVertexOnCurve(edges[0], 0.25)
    vertex2 = geompy.MakeVertexOnCurve(edges[0], 0.75)
    basicCoords1 = geompy.PointCoordinates(vertex1)
    basicCoords2 = geompy.PointCoordinates(vertex2)
    tgtEdges = geompy.GetInPlace(Partition_2, lEdges)
    edges = geompy.SubShapeAll(tgtEdges, geompy.ShapeType["EDGE"])
    alongx_raff3 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
    alongy_raff3 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])

    if (abs(basicCoords2[0] - basicCoords1[0])) > 5:
        # Selon X
        edges_x_raff3 = edges_x_raff3 + edges
    if (abs(basicCoords2[1] - basicCoords1[1])) > 5:
        # Selon Y
        edges_y_raff3 = edges_y_raff3 + edges
geompy.UnionList(alongx_raff3, edges_x_raff3)
geompy.UnionList(alongy_raff3, edges_y_raff3)

# Récuperation des edges du domaine
edge_dom = []
list_edges_dom = geompy.GetShapesOnShape(
    domaine, Partition_2, geompy.ShapeType["EDGE"], GEOM.ST_ONIN
)
for ledge in list_edges_dom:
    g_edge_dom = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
    edge_dom.append(ledge)
geompy.UnionList(g_edge_dom, edge_dom)

# Récuperation des zones de raffinement parc raff4
edges_x_raff4: List = []
edges_y_raff4: List = []
lListEdges = geompy.Propagate(Union_raff4)
for lEdges in lListEdges:
    edges = geompy.SubShapeAll(lEdges, geompy.ShapeType["EDGE"])
    vertex1 = geompy.MakeVertexOnCurve(edges[0], 0.25)
    vertex2 = geompy.MakeVertexOnCurve(edges[0], 0.75)
    basicCoords1 = geompy.PointCoordinates(vertex1)
    basicCoords2 = geompy.PointCoordinates(vertex2)
    tgtEdges = geompy.GetInPlace(Partition_2, lEdges)
    edges = geompy.SubShapeAll(tgtEdges, geompy.ShapeType["EDGE"])
    alongx_raff4 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
    alongy_raff4 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])

    if (abs(basicCoords2[0] - basicCoords1[0])) > 5:
        # Selon X
        edges_x_raff4 = edges_x_raff4 + edges
    if (abs(basicCoords2[1] - basicCoords1[1])) > 5:
        # Selon Y
        edges_y_raff4 = edges_y_raff4 + edges
geompy.UnionList(alongx_raff4, edges_x_raff4)
geompy.UnionList(alongy_raff4, edges_y_raff4)

# Récuperation des zones de raffinement parc raff5
edges_x_raff5: List = []
edges_y_raff5: List = []
lListEdges = geompy.Propagate(Union_raff5)
for lEdges in lListEdges:
    edges = geompy.SubShapeAll(lEdges, geompy.ShapeType["EDGE"])
    vertex1 = geompy.MakeVertexOnCurve(edges[0], 0.25)
    vertex2 = geompy.MakeVertexOnCurve(edges[0], 0.75)
    basicCoords1 = geompy.PointCoordinates(vertex1)
    basicCoords2 = geompy.PointCoordinates(vertex2)
    tgtEdges = geompy.GetInPlace(Partition_2, lEdges)
    edges = geompy.SubShapeAll(tgtEdges, geompy.ShapeType["EDGE"])
    alongx_raff5 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
    alongy_raff5 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])

    if (abs(basicCoords2[0] - basicCoords1[0])) > 5:
        # Selon X
        edges_x_raff5 = edges_x_raff5 + edges
    if (abs(basicCoords2[1] - basicCoords1[1])) > 5:
        # Selon Y
        edges_y_raff5 = edges_y_raff5 + edges
geompy.UnionList(alongx_raff5, edges_x_raff5)
geompy.UnionList(alongy_raff5, edges_y_raff5)

# Récuperation des zones de raffinement parc raff6
edges_x_raff6: List = []
edges_y_raff6: List = []
lListEdges = geompy.Propagate(Union_raff6)
for lEdges in lListEdges:
    edges = geompy.SubShapeAll(lEdges, geompy.ShapeType["EDGE"])
    vertex1 = geompy.MakeVertexOnCurve(edges[0], 0.25)
    vertex2 = geompy.MakeVertexOnCurve(edges[0], 0.75)
    basicCoords1 = geompy.PointCoordinates(vertex1)
    basicCoords2 = geompy.PointCoordinates(vertex2)
    tgtEdges = geompy.GetInPlace(Partition_2, lEdges)
    edges = geompy.SubShapeAll(tgtEdges, geompy.ShapeType["EDGE"])
    alongx_raff6 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
    alongy_raff6 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])

    if (abs(basicCoords2[0] - basicCoords1[0])) > 5:
        # Selon X
        edges_x_raff6 = edges_x_raff6 + edges
    if (abs(basicCoords2[1] - basicCoords1[1])) > 5:
        # Selon Y
        edges_y_raff6 = edges_y_raff6 + edges
geompy.UnionList(alongx_raff6, edges_x_raff6)
geompy.UnionList(alongy_raff6, edges_y_raff6)

# Récuperation des zones de raffinement parc raff7
edges_x_raff7: List = []
edges_y_raff7: List = []
lListEdges = geompy.Propagate(Union_raff7)
for lEdges in lListEdges:
    edges = geompy.SubShapeAll(lEdges, geompy.ShapeType["EDGE"])
    vertex1 = geompy.MakeVertexOnCurve(edges[0], 0.25)
    vertex2 = geompy.MakeVertexOnCurve(edges[0], 0.75)
    basicCoords1 = geompy.PointCoordinates(vertex1)
    basicCoords2 = geompy.PointCoordinates(vertex2)
    tgtEdges = geompy.GetInPlace(Partition_2, lEdges)
    edges = geompy.SubShapeAll(tgtEdges, geompy.ShapeType["EDGE"])
    alongx_raff7 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
    alongy_raff7 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])

    if (abs(basicCoords2[0] - basicCoords1[0])) > 5:
        # Selon X
        edges_x_raff7 = edges_x_raff7 + edges
    if (abs(basicCoords2[1] - basicCoords1[1])) > 5:
        # Selon Y
        edges_y_raff7 = edges_y_raff7 + edges
geompy.UnionList(alongx_raff7, edges_x_raff7)
geompy.UnionList(alongy_raff7, edges_y_raff7)

# pour extrusion du maillage
marge = 10
end_point = geompy.MakeVertex(0, 0, extrusion1_max_height + marge)
edge_wire = geompy.MakeEdge(origin, end_point)
wire = geompy.MakeWire([edge_wire], 1e-07)


if args.damping_layer > 0:
    end_point2 = geompy.MakeVertex(0, 0, extrusion2_height)
    end_point3 = geompy.MakeVertex(0, 0, Lz)
    edge_wire3 = geompy.MakeEdge(end_point2, end_point3)
    wire3 = geompy.MakeWire([edge_wire3], 1e-07)
else:
    end_point2 = geompy.MakeVertex(0, 0, Lz)
#
edge_wire2 = geompy.MakeEdge(end_point, end_point2)
wire2 = geompy.MakeWire([edge_wire2], 1e-07)

Point = geompy.MakeVertex(-Lx / 2, 0, 0)
Plane_x0 = geompy.MakePlane(Point, OX, 10000)

Point = geompy.MakeVertex(Lx / 2, 0, 0)
Plane_x1 = geompy.MakePlane(Point, OX, 10000)

Point = geompy.MakeVertex(0, -Ly / 2, 0)
Plane_y0 = geompy.MakePlane(Point, OY, 10000)

Point = geompy.MakeVertex(0, Ly / 2, 0)
Plane_y1 = geompy.MakePlane(Point, OY, 10000)

edge_bande = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(edge_bande, [9, 13])
edge_cercle_int = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(edge_cercle_int, [11, 16])
edge_cercle_ext = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(edge_cercle_ext, [4, 19, 17, 7])
face_bande = geompy.CreateGroup(Partition_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(face_bande, [2, 14])

ppts = []
ppts2 = {}
for i in range(nombre_turbines):
    ppts.append(np.array(xy_turbines[:, i]))
    ppts2[i] = list([np.array(xy_turbines[:, i])])
ppts = np.array(ppts)

# Stocker coordonnees dans fichier pour les simus
# xm,ym=ppts.T
# ppts=np.array(ppts)
# with open("turbines","w") as fic:
#  fic.write("%i"%len(xm)+os.linesep)
#  for i,(x,y) in enumerate(zip(xm,ym)):
#    fic.write("%i"%(i+1,)+os.linesep)
#    fic.write("%.2f"%x+os.linesep)
#    fic.write("%.2f"%y+os.linesep)
#    fic.write("%.2f"%hm+os.linesep)

###
### SMESH component
###

import SALOMEDS
import SMESH
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

Maillage_1 = smesh.Mesh(Partition_2, "Maillage_1")

# zone raff1
Regular_1D_raff1_x = Maillage_1.Segment(geom=alongx_raff1)
Number_of_Segments_raff1_x = Regular_1D_raff1_x.LocalLength(tm, None, 1e-07)

Regular_1D_raff1_y = Maillage_1.Segment(geom=alongy_raff1)
Number_of_Segments_raff1_y = Regular_1D_raff1_y.LocalLength(tm, None, 1e-07)

Regular_1D_raff1_x = Regular_1D_raff1_x.GetSubMesh()
Regular_1D_raff1_y = Regular_1D_raff1_y.GetSubMesh()

# zone raff2
Regular_1D_raff2_x = Maillage_1.Segment(geom=alongx_raff2)
Number_of_Segments_raff2_x = Regular_1D_raff2_x.LocalLength(traff2, None, 1e-07)

Regular_1D_raff2_y = Maillage_1.Segment(geom=alongy_raff2)
Number_of_Segments_raff2_y = Regular_1D_raff2_y.LocalLength(traff2, None, 1e-07)

Regular_1D_raff2_x = Regular_1D_raff2_x.GetSubMesh()
Regular_1D_raff2_y = Regular_1D_raff2_y.GetSubMesh()

# zone raff3
Regular_1D_raff3_x = Maillage_1.Segment(geom=alongx_raff3)
Number_of_Segments_raff3_x = Regular_1D_raff3_x.LocalLength(traff3, None, 1e-07)

Regular_1D_raff3_y = Maillage_1.Segment(geom=alongy_raff3)
Number_of_Segments_raff3_y = Regular_1D_raff3_y.LocalLength(traff3, None, 1e-07)

Regular_1D_raff3_x = Regular_1D_raff3_x.GetSubMesh()
Regular_1D_raff3_y = Regular_1D_raff3_y.GetSubMesh()

# zone raff4
Regular_1D_raff4_x = Maillage_1.Segment(geom=alongx_raff4)
Number_of_Segments_raff4_x = Regular_1D_raff4_x.LocalLength(traff4, None, 1e-07)

Regular_1D_raff4_y = Maillage_1.Segment(geom=alongy_raff4)
Number_of_Segments_raff4_y = Regular_1D_raff4_y.LocalLength(traff4, None, 1e-07)

Regular_1D_raff4_x = Regular_1D_raff4_x.GetSubMesh()
Regular_1D_raff4_y = Regular_1D_raff4_y.GetSubMesh()

# zone raff5
Regular_1D_raff5_x = Maillage_1.Segment(geom=alongx_raff5)
Number_of_Segments_raff5_x = Regular_1D_raff5_x.LocalLength(traff5, None, 1e-07)

Regular_1D_raff5_y = Maillage_1.Segment(geom=alongy_raff5)
Number_of_Segments_raff5_y = Regular_1D_raff5_y.LocalLength(traff5, None, 1e-07)

Regular_1D_raff5_x = Regular_1D_raff5_x.GetSubMesh()
Regular_1D_raff5_y = Regular_1D_raff5_y.GetSubMesh()

# zone raff6
Regular_1D_raff6_x = Maillage_1.Segment(geom=alongx_raff6)
Number_of_Segments_raff6_x = Regular_1D_raff6_x.LocalLength(traff6, None, 1e-07)

Regular_1D_raff6_y = Maillage_1.Segment(geom=alongy_raff6)
Number_of_Segments_raff6_y = Regular_1D_raff6_y.LocalLength(traff6, None, 1e-07)

Regular_1D_raff6_x = Regular_1D_raff6_x.GetSubMesh()
Regular_1D_raff6_y = Regular_1D_raff6_y.GetSubMesh()

# zone raff7
Regular_1D_raff7_x = Maillage_1.Segment(geom=alongx_raff7)
Number_of_Segments_raff7_x = Regular_1D_raff7_x.LocalLength(tc, None, 1e-07)

Regular_1D_raff7_y = Maillage_1.Segment(geom=alongy_raff7)
Number_of_Segments_raff7_y = Regular_1D_raff7_y.LocalLength(tc, None, 1e-07)

Regular_1D_raff7_x = Regular_1D_raff7_x.GetSubMesh()
Regular_1D_raff7_y = Regular_1D_raff7_y.GetSubMesh()

paroi = Maillage_1.GroupOnGeom(g_edge_dom, "edgex0_dom", SMESH.EDGE)

if dom == 0:
    # nb = 2*pi*(max(Lx/2,Ly/2)-lb)/(4*tc)
    nb = 2 * pi * (max(Lx / 2, Ly / 2)) / (4 * tc)
    nb2 = lb / tc
    Regular_1D_bande = Maillage_1.Segment(geom=edge_bande)
    Number_of_Segments_bande = Regular_1D_bande.NumberOfSegments(ceil(nb2))
    Regular_1D_cercleint = Maillage_1.Segment(geom=edge_cercle_int)
    Number_of_Segments_cercleint = Regular_1D_cercleint.NumberOfSegments(2 * ceil(nb))
    Regular_1D_cercleext = Maillage_1.Segment(geom=edge_cercle_ext)
    Number_of_Segments_cercleext = Regular_1D_cercleext.NumberOfSegments(ceil(nb))
    Quadrangle_2D_1 = Maillage_1.Quadrangle(
        algo=smeshBuilder.QUADRANGLE, geom=face_bande
    )

tmax = max(tc, tv2)
tmin = min(tm, tv1)
GMSH_2D = Maillage_1.Triangle(algo=smeshBuilder.GMSH_2D)
Gmsh_Parameters = GMSH_2D.Parameters()
Gmsh_Parameters.Set2DAlgo(4)
Gmsh_Parameters.SetRecombineAll(1)
Gmsh_Parameters.SetMinSize(tmin)
Gmsh_Parameters.SetMaxSize(tmax)
# Gmsh_Parameters.SetSmouthSteps( 100 )
Gmsh_Parameters.SetSubdivAlgo(0)
Gmsh_Parameters.SetRemeshAlgo(0)
Gmsh_Parameters.SetMeshCurvatureSize(0)
Gmsh_Parameters.SetRecomb2DAlgo(1)
Gmsh_Parameters.SetIs2d(1)
isDone = Maillage_1.Compute()

# path pour l'extrusion

path_meshz = smesh.Mesh(wire)
Regular_1D_z = path_meshz.Segment()
Local_Length_z = Regular_1D_z.LocalLength(tv1)
# Number_of_Segments_z = Regular_1D_z.StartEndLength(7,100)
# Number_of_Segments_z = Regular_1D_z.GeometricProgression(5,1.05,[])
# Local_Length_1.SetPrecision( 1e-07 )
isDone = path_meshz.Compute()
smesh.SetName(path_meshz.GetMesh(), "path_1")

path_meshz2 = smesh.Mesh(wire2)
Regular_1D_z2 = path_meshz2.Segment()
# Local_Length_z2 =  Regular_1D_z2.GeometricProgression(6,1.15,[])
Number_of_Segments_z2 = Regular_1D_z2.StartEndLength(tv1, tv2)
isDone = path_meshz2.Compute()
smesh.SetName(path_meshz2.GetMesh(), "path_2")

if args.damping_layer > 0.0:
    path_meshz3 = smesh.Mesh(wire3)
    Regular_1D_z3 = path_meshz3.Segment()
    # Local_Length_z2 =  Regular_1D_z2.GeometricProgression(6,1.15,[])
    Number_of_Segments_z3 = Regular_1D_z3.LocalLength(tv2)
    isDone = path_meshz3.Compute()
    smesh.SetName(path_meshz3.GetMesh(), "path_3")
    path = smesh.Concatenate(
        [path_meshz.GetMesh(), path_meshz2.GetMesh(), path_meshz3.GetMesh()],
        1,
        1,
        1e-05,
        False,
    )
else:
    path = smesh.Concatenate(
        [path_meshz.GetMesh(), path_meshz2.GetMesh()], 1, 1, 1e-05, False
    )

smesh.SetName(path.GetMesh(), "path")

# extrusion
coords = geompy.PointCoordinates(origin)
start_node = path_meshz.FindNodeClosestTo(coords[0], coords[1], coords[2])
isDone = Maillage_1.Compute()
(Group, error) = Maillage_1.ExtrusionAlongPathObjects(
    [Maillage_1], [], [Maillage_1], path, None, 1, 0, [], 0, 0, [0, 0, 0], 1, [], 0
)
smesh.SetName(Maillage_1.GetMesh(), "Maillage_1")

listGroup = []
for i in range(len(Group)):
    if Group[i].GetType() == SMESH.FACE:
        listGroup.append(Group[i])

Group[0].SetName("Bords")

Sol = Maillage_1.GroupOnGeom(Partition_2, "Partition_2", SMESH.FACE)
Sol.SetName("Sol")

# groupe de paroi : toutes les faces externes sauf celles déjà dans des groupes
faces_externes = Maillage_1.CreateEmptyGroup(SMESH.FACE, "faces_externes")
nbAdd = faces_externes.AddFrom(Maillage_1.GetMesh())
listGroup.append(Sol)
Sommet = Maillage_1.GetMesh().CutListOfGroups([faces_externes], listGroup, "Sommet")
Maillage_1.RemoveGroup(faces_externes)

Maillage_1.ExportMED(medFile, 0)

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()
