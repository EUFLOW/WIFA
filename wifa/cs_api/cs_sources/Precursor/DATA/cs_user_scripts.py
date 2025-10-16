#!/usr/bin/env python

# -------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

# -------------------------------------------------------------------------------

import os

# ===============================================================================
# Local functions
# ===============================================================================

# -------------------------------------------------------------------------------
# Process the command line
# -------------------------------------------------------------------------------


def process_command_line(argv):
    """
    Processes the passed command line arguments.
    """
    from optparse import OptionParser

    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option(
        "--meteo",
        "--meteo-file",
        dest="MeteoFile",
        type="string",
        help="Name of meteo file",
    )

    (options, args) = parser.parse_args(argv)

    return options


# ===============================================================================
# Local functions
# ===============================================================================

# ===============================================================================
# Defining parameters for a calculation domain
# ===============================================================================


def domain_prepare_data_add(domain):
    """
    Additional steps to prepare data
    (called in data preparation stage, between copy of files
    in DATA and copy of link of restart files as defined by domain).
    """

    return


# -------------------------------------------------------------------------------


def domain_copy_results_add(domain):
    """
    Additional steps to copy results or cleanup execution directory
    (called at beginning of data copy stage).
    """

    return


# -------------------------------------------------------------------------------


def define_domain_parameters(domain):
    """
    Define domain execution parameters.
    """

    # Path for cronos librairies
    # Don't forget to change if used on gaia or another cluster
    domain.compile_cflags = "-I/software/rd/saturne/usr/include/gdal"
    domain.compile_cxxflags = None
    domain.compile_fcflags = None
    domain.compile_libs = "-lgdal -lproj"

    if domain.kw_args:
        options = process_command_line(domain.kw_args)

        from code_saturne.model.AtmosphericFlowsModel import AtmosphericFlowsModel
        from code_saturne.model.XMLengine import Case
        from code_saturne.model.XMLinitialize import XMLinit

        fp = os.path.join(domain.exec_dir, "setup.xml")
        case = Case(package=domain.package, file_name=fp)
        case["xmlfile"] = fp
        case.xmlCleanAllBlank(case.xmlRootNode())
        XMLinit(case).initialize()

        if options.MeteoFile:
            t = AtmosphericFlowsModel(case)
            fullPathMeteo = options.MeteoFile
            t.setMeteoDataFileName(fullPathMeteo)

        case.xmlSaveDocument()

    return


# -------------------------------------------------------------------------------
# End
# -------------------------------------------------------------------------------
