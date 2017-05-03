# Copyright (C) 2015-2017: The University of Edinburgh
#                 Authors: Craig Warren and Antonis Giannopoulos
#
# This file is part of gprMax.
#
# gprMax is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gprMax is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gprMax.  If not, see <http://www.gnu.org/licenses/>.

import datetime
import itertools
import os
import psutil
import sys

from colorama import init, Fore, Style

import numpy as np
from terminaltables import AsciiTable
from tqdm import tqdm

from .exceptions import GeneralError

from .fields_outputs import write_hdf5_outputfile

from .grid import FDTDGrid
from .grid import dispersion_analysis
from .input_cmds_geometry import process_geometrycmds
from .input_cmds_file import process_python_include_code
from .input_cmds_file import write_processed_file
from .input_cmds_file import check_cmd_names
from .input_cmds_multiuse import process_multicmds
from .input_cmds_singleuse import process_singlecmds
from .materials import Material, process_materials
from .pml import build_pmls
from .utilities import get_terminal_width
from .utilities import human_size
from .yee_cell_build import build_electric_components
from .yee_cell_build import build_magnetic_components
from .solvers import CPUSolver

init()


def run_model(args, currentmodelrun, modelend, numbermodelruns, inputfile,
              usernamespace):
    """Runs a model - processes the input file; builds the Yee cells;
    calculates update coefficients; runs main FDTD loop.

    Args:
        args (dict): Namespace with command line arguments
        currentmodelrun (int): Current model run number.
        modelend (int): Number of last model to run.
        numbermodelruns (int): Total number of model runs.
        inputfile (object): File object for the input file.
        usernamespace (dict): Namespace that can be accessed by user
                in any Python code blocks in input file.

    Returns:
        tsolve (int): Length of time (seconds) of main FDTD calculations
    """

    # Monitor memory usage
    p = psutil.Process()

    # Declare variable to hold FDTDGrid class
    global G

    # Used for naming geometry and output files
    appendmodelnumber = '' if numbermodelruns == 1 and not args.task and not args.restart else str(currentmodelrun)

    # Normal model reading/building process; bypassed if geometry information to be reused
    if 'G' not in globals():

        # Initialise an instance of the FDTDGrid class
        G = FDTDGrid()

        G.inputfilename = os.path.split(inputfile.name)[1]
        G.inputdirectory = os.path.dirname(os.path.abspath(inputfile.name))
        inputfilestr = '\n--- Model {}/{}, input file: {}'.format(currentmodelrun, modelend, inputfile.name)
        print(Fore.GREEN + '{} {}\n'.format(inputfilestr, '-' * (get_terminal_width() - 1 - len(inputfilestr))) + Style.RESET_ALL)

        # Add the current model run to namespace that can be accessed by
        # user in any Python code blocks in input file
        usernamespace['current_model_run'] = currentmodelrun

        # Read input file and process any Python and include file commands
        processedlines = process_python_include_code(inputfile, usernamespace)

        # Print constants/variables in user-accessable namespace
        uservars = ''
        for key, value in sorted(usernamespace.items()):
            if key != '__builtins__':
                uservars += '{}: {}, '.format(key, value)
        print('Constants/variables used/available for Python scripting: {{{}}}\n'.format(uservars[:-2]))

        # Write a file containing the input commands after Python or include file commands have been processed
        if args.write_processed:
            write_processed_file(processedlines, appendmodelnumber, G)

        # Check validity of command names and that essential commands are present
        singlecmds, multicmds, geometry = check_cmd_names(processedlines)

        # Create built-in materials
        m = Material(0, 'pec')
        m.se = float('inf')
        m.type = 'builtin'
        m.averagable = False
        G.materials.append(m)
        m = Material(1, 'free_space')
        m.type = 'builtin'
        G.materials.append(m)

        # Process parameters for commands that can only occur once in the model
        process_singlecmds(singlecmds, G)

        # Process parameters for commands that can occur multiple times in the model
        print()
        process_multicmds(multicmds, G)

        # Initialise an array for volumetric material IDs (solid), boolean
        # arrays for specifying materials not to be averaged (rigid),
        # an array for cell edge IDs (ID)
        G.initialise_geometry_arrays()

        # Initialise arrays for the field components
        G.initialise_field_arrays()

        # Process geometry commands in the order they were given
        process_geometrycmds(geometry, G)

        # Build the PMLs and calculate initial coefficients
        print()
        if all(value == 0 for value in G.pmlthickness.values()):
            if G.messages:
                print('PML boundaries: switched off')
            pass  # If all the PMLs are switched off don't need to build anything
        else:
            if G.messages:
                if all(value == G.pmlthickness['x0'] for value in G.pmlthickness.values()):
                    pmlinfo = str(G.pmlthickness['x0']) + ' cells'
                else:
                    pmlinfo = ''
                    for key, value in G.pmlthickness.items():
                        pmlinfo += '{}: {} cells, '.format(key, value)
                    pmlinfo = pmlinfo[:-2]
                print('PML boundaries: {}'.format(pmlinfo))
            pbar = tqdm(total=sum(1 for value in G.pmlthickness.values() if value > 0), desc='Building PML boundaries', ncols=get_terminal_width() - 1, file=sys.stdout, disable=G.tqdmdisable)
            build_pmls(G, pbar)
            pbar.close()

        # Build the model, i.e. set the material properties (ID) for every edge
        # of every Yee cell
        print()
        pbar = tqdm(total=2, desc='Building main grid', ncols=get_terminal_width() - 1, file=sys.stdout, disable=G.tqdmdisable)
        build_electric_components(G.solid, G.rigidE, G.ID, G)
        pbar.update()
        build_magnetic_components(G.solid, G.rigidH, G.ID, G)
        pbar.update()
        pbar.close()

        # Process any voltage sources (that have resistance) to create a new
        # material at the source location
        for voltagesource in G.voltagesources:
            voltagesource.create_material(G)

        # Initialise arrays of update coefficients to pass to update functions
        G.initialise_std_update_coeff_arrays()

        # Initialise arrays of update coefficients and temporary values if
        # there are any dispersive materials
        if Material.maxpoles != 0:
            G.initialise_dispersive_arrays()

        # Process complete list of materials - calculate update coefficients,
        # store in arrays, and build text list of materials/properties
        materialsdata = process_materials(G)
        if G.messages:
            print('\nMaterials:')
            materialstable = AsciiTable(materialsdata)
            materialstable.outer_border = False
            materialstable.justify_columns[0] = 'right'
            print(materialstable.table)

        # Check to see if numerical dispersion might be a problem
        results = dispersion_analysis(G)
        if not results['waveform']:
            print(Fore.RED + "\nWARNING: Numerical dispersion analysis not carried out as either no waveform detected or waveform does not fit within specified time window and is therefore being truncated." + Style.RESET_ALL)
        elif results['N'] < G.mingridsampling:
            raise GeneralError("Non-physical wave propagation: Material '{}' has wavelength sampled by {} cells, less than required minimum for physical wave propagation. Maximum significant frequency estimated as {:g}Hz".format(results['material'].ID, results['N'], results['maxfreq']))
        elif results['deltavp'] and np.abs(results['deltavp']) > G.maxnumericaldisp:
            print(Fore.RED + "\nWARNING: Potentially significant numerical dispersion. Estimated largest physical phase-velocity error is {:.2f}% in material '{}' whose wavelength sampled by {} cells. Maximum significant frequency estimated as {:g}Hz".format(results['deltavp'], results['material'].ID, results['N'], results['maxfreq']) + Style.RESET_ALL)
        elif results['deltavp'] and G.messages:
            print("\nNumerical dispersion analysis: estimated largest physical phase-velocity error is {:.2f}% in material '{}' whose wavelength sampled by {} cells. Maximum significant frequency estimated as {:g}Hz".format(results['deltavp'], results['material'].ID, results['N'], results['maxfreq']))

    # If geometry information to be reused between model runs
    else:
        inputfilestr = '\n--- Model {}/{}, input file (not re-processed, i.e. geometry fixed): {}'.format(currentmodelrun, modelend, inputfile.name)
        print(Fore.GREEN + '{} {}\n'.format(inputfilestr, '-' * (get_terminal_width() - 1 - len(inputfilestr))) + Style.RESET_ALL)

        # Clear arrays for field components
        G.initialise_field_arrays()

        # Clear arrays for fields in PML
        for pml in G.pmls:
            pml.initialise_field_arrays()

    # Adjust position of simple sources and receivers if required
    if G.srcsteps[0] != 0 or G.srcsteps[1] != 0 or G.srcsteps[2] != 0:
        for source in itertools.chain(G.hertziandipoles, G.magneticdipoles):
            if currentmodelrun == 1:
                if source.xcoord + G.srcsteps[0] * modelend < 0 or source.xcoord + G.srcsteps[0] * modelend > G.nx or source.ycoord + G.srcsteps[1] * modelend < 0 or source.ycoord + G.srcsteps[1] * modelend > G.ny or source.zcoord + G.srcsteps[2] * modelend < 0 or source.zcoord + G.srcsteps[2] * modelend > G.nz:
                    raise GeneralError('Source(s) will be stepped to a position outside the domain.')
            source.xcoord = source.xcoordorigin + (currentmodelrun - 1) * G.srcsteps[0]
            source.ycoord = source.ycoordorigin + (currentmodelrun - 1) * G.srcsteps[1]
            source.zcoord = source.zcoordorigin + (currentmodelrun - 1) * G.srcsteps[2]
    if G.rxsteps[0] != 0 or G.rxsteps[1] != 0 or G.rxsteps[2] != 0:
        for receiver in G.rxs:
            if currentmodelrun == 1:
                if receiver.xcoord + G.rxsteps[0] * modelend < 0 or receiver.xcoord + G.rxsteps[0] * modelend > G.nx or receiver.ycoord + G.rxsteps[1] * modelend < 0 or receiver.ycoord + G.rxsteps[1] * modelend > G.ny or receiver.zcoord + G.rxsteps[2] * modelend < 0 or receiver.zcoord + G.rxsteps[2] * modelend > G.nz:
                    raise GeneralError('Receiver(s) will be stepped to a position outside the domain.')
            receiver.xcoord = receiver.xcoordorigin + (currentmodelrun - 1) * G.rxsteps[0]
            receiver.ycoord = receiver.ycoordorigin + (currentmodelrun - 1) * G.rxsteps[1]
            receiver.zcoord = receiver.zcoordorigin + (currentmodelrun - 1) * G.rxsteps[2]

    # Write files for any geometry views and geometry object outputs
    if not (G.geometryviews or G.geometryobjectswrite) and args.geometry_only:
        print(Fore.RED + '\nWARNING: No geometry views or geometry objects to output found.' + Style.RESET_ALL)
    if G.geometryviews:
        print()
        for i, geometryview in enumerate(G.geometryviews):
            geometryview.set_filename(appendmodelnumber, G)
            pbar = tqdm(total=geometryview.datawritesize, unit='byte', unit_scale=True, desc='Writing geometry view file {}/{}, {}'.format(i + 1, len(G.geometryviews), os.path.split(geometryview.filename)[1]), ncols=get_terminal_width() - 1, file=sys.stdout, disable=G.tqdmdisable)
            geometryview.write_vtk(G, pbar)
            pbar.close()
    if G.geometryobjectswrite:
        for i, geometryobject in enumerate(G.geometryobjectswrite):
            pbar = tqdm(total=geometryobject.datawritesize, unit='byte', unit_scale=True, desc='Writing geometry object file {}/{}, {}'.format(i + 1, len(G.geometryobjectswrite), os.path.split(geometryobject.filename)[1]), ncols=get_terminal_width() - 1, file=sys.stdout, disable=G.tqdmdisable)
            geometryobject.write_hdf5(G, pbar)
            pbar.close()

    # If only writing geometry information
    if args.geometry_only:
        tsolve = 0

    # Run simulation
    else:
        # Prepare any snapshot files
        for snapshot in G.snapshots:
            snapshot.prepare_vtk_imagedata(appendmodelnumber, G)

        # Output filename
        inputfileparts = os.path.splitext(os.path.join(G.inputdirectory, G.inputfilename))
        outputfile = inputfileparts[0] + appendmodelnumber + '.out'
        print('\nOutput file: {}\n'.format(outputfile))

        # Main FDTD solving functions for either CPU or GPU
        desc = 'Running simulation, model {}/{}'.format(str(currentmodelrun),
                                                        str(modelend))
        iterations = tqdm(
                          range(G.iterations),
                          desc=desc,
                          ncols=get_terminal_width() - 1,
                          file=sys.stdout,
                          disable=G.tqdmdisable)

        cpusolver = CPUSolver(G, iterations)
        tsolve = cpusolver.solve()

        # Write an output file in HDF5 format
        write_hdf5_outputfile(outputfile, G.Ex, G.Ey, G.Ez, G.Hx, G.Hy, G.Hz, G)

        if G.messages:
            print('Memory (RAM) used: ~{}'.format(human_size(p.memory_info().rss)))
            print('Solving time [HH:MM:SS]: {}'.format(datetime.timedelta(seconds=tsolve)))

        # If geometry information to be reused between model runs then FDTDGrid
        # class instance must be global so that it persists
        if not args.geometry_fixed:
            del G

    return tsolve
