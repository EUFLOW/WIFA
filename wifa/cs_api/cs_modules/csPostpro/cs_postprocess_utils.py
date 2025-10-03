#!/usr/bin/env python3
import sys

import matplotlib.tri as tri
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy as vtk_to_np


def load_ensight_data(path, case):
    """
    load the ensight file
    """
    ens = vtk.vtkGenericEnSightReader()
    ens.SetFilePath(path)
    ens.SetCaseFileName(case)
    ens.ReadAllVariablesOn()
    ens.Update()
    return ens


def read_ensight_data_at_time(ens, inst=None):
    """
    read the results for a given time inst
    """
    try:
        times = ens.GetTimeSets().GetItem(0)
        if inst > times.GetSize() - 1 or inst < 0 or inst is None:
            print(
                "Time Exceeds max or is negative. Read time %i" % (times.GetSize() - 1)
            )
            inst = times.GetSize() - 1
        ens.SetTimeValue(times.GetTuple1(inst))
    except AttributeError:
        print("Error reading the required time step")
        pass
    ens.Update()
    data = ens.GetOutput()
    if data.GetClassName() == "vtkMultiBlockDataSet":
        data = data.GetBlock(0)

    return data


def ensight_time_list(ens):
    """
    read the time list
    """
    try:
        times = ens.GetTimeSets().GetItem(0)
        return list(range(times.GetNumberOfTuples())), [
            times.GetTuple1(i) for i in range(times.GetNumberOfTuples())
        ]
    except AttributeError:
        return []


def ensight_var_list(data):
    """
    read the list of variable names
    """
    fields = getattr(data, "GetCellData")()
    varname_list = [fields.GetArrayName(i) for i in range(fields.GetNumberOfArrays())]
    return varname_list


def get_field_from_ensight(data, varname):
    """
    extract the field 'varname' from ensight 'data'
    """
    field = getattr(data, "GetCellData")().GetArray(varname)
    if field is not None:
        tab = vtk_to_np(field)
        tab = tab.copy()
        return tab
    else:
        return None


def get_point_coords_from_ensight(data):
    """
    get the coordinates of mesh points
    """
    return np.array(
        [data.GetPoint(i) for i in range(data.GetNumberOfPoints())]
    ).transpose()


def extract_plane_from_ensight(data, origin=(0.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)):
    """
    extract plane of center 'origin' orthogonally to 'normal'
    """
    plane = vtk.vtkPlane()
    plane.SetOrigin(*origin)
    plane.SetNormal(*normal)
    extract = vtk.vtkCutter()
    extract.SetInputData(data)
    extract.SetCutFunction(plane)
    extract.Update()
    plane_cut = extract.GetOutput()

    if not plane_cut.GetNumberOfCells():
        print("Empty plane - PASS")
        return None

    plane_cut.GetCellData().SetActiveScalars("Velocity")

    return plane_cut


def extract_saturne_triangulation(slicexy, normal):
    """extract the triangulation of a horizontal slice of data"""
    triangles = slicexy.GetPolys().GetData()
    npts = slicexy.GetPoints().GetNumberOfPoints()
    ntri = int(triangles.GetNumberOfTuples() / 4)
    x, y, z = get_point_coords_from_ensight(slicexy)  # points coordinates

    if normal == (0, 0, 1):
        x1 = x
        y1 = y
    if normal == (1, 0, 0):
        x1 = y
        y1 = z
    if normal == (0, 1, 0):
        x1 = x
        y1 = z

    # Get the triangulation of the mesh
    triang = np.zeros((ntri, 3))
    for i in range(0, ntri):
        triang[i, 0] = triangles.GetTuple(4 * i + 1)[0]
        triang[i, 1] = triangles.GetTuple(4 * i + 2)[0]
        triang[i, 2] = triangles.GetTuple(4 * i + 3)[0]

    triangulation = tri.Triangulation(x1, y1, triang)

    return x1, y1, triang, triangulation


def get_cell_centers(data):
    """
    Get the coordinates of the cell centers
    """
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputData(data)
    cell_centers.Update()
    centers = cell_centers.GetOutput()
    points = centers.GetPoints()

    return points
