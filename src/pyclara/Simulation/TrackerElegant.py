import h5py

import pyclara
import sdds
import subprocess
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import os
import re


class Elegant_runner():
    """
    Runs Elegant beam tracking simulations as part of an S2E pipeline.

    Handles reading input files, converting particle distributions,
    writing Elegant input, running the simulation, and returning output
    particles in OpenPMD HDF5 format.

    Parameters
    ----------
    inputdir : str
        Path to the directory containing Elegant input files (.lte, .ele).
    name : str
        Name of the simulation and lattice (e.g. 'FEBE'). Used to locate
        files named {name}.lte and {name}.ele.
    run_dir : str
        Path to the directory where simulation files will be written and run.

    Attributes
    ----------
    name : str
        Simulation name.
    run_dir : str
        Directory where simulation is executed.
    elegant_lte : str
        Contents of the lattice file ({name}.lte).
    elegant_ele : str
        Contents of the Elegant run file ({name}.ele).
    wakefiles : dict
        Dictionary of {filename: content} for any wake files referenced
        in the lattice.
        
    Examples
    --------
    # e = Elegant_runner('/path/to/elegant_input', 'FEBE', '/path/to/output')
    # e.set_input_h5(h5py.File('input.h5', 'r'))
    # e.write_elegant_input()
    # e.run()
    # output = e.get_output()
    """
    def __init__(self, inputdir, name, run_dir):
        self.name = name
        self.run_dir = run_dir
        elegant_lte = f'{inputdir}/{name}.lte'
        elegant_ele = f'{inputdir}/{name}.ele'
        self.ele_path = elegant_ele
        with open(elegant_lte, 'r') as r:
            self.elegant_lte = r.read()
        with open(elegant_ele, 'r') as r:
            self.elegant_ele = r.read()
        self.wakefiles = {}
        seen = set()
        l = pyclara.Converters.elegant_lte_loader(elegant_lte)
        for k in l:
            if "WAKEFILE" in l[k]:
                n = l[k]["WAKEFILE"]
                if n not in seen:
                    file = f'{inputdir}/{n}'
                    with open(file, 'r') as f:
                        self.wakefiles[n] = f.read()
        self.cmd = [
            "/lib64/mpich/bin/mpirun",
            "-np", "6",
            "/elegant/elegant/bin/Linux-x86_64/Pelegant",
            f"{self.name}.ele"
        ]

    def set_input_sdds(self, input_particle):
        """
                Set the input particle distribution from an SDDS file.

                Parameters
                ----------
                input_particle : str
                    Path to the input SDDS file.
        """
        self.sdds_input_path = input_particle
        self.sdds_input = sdds.load(self.sdds_input_path)
        self.detect = False

    def set_input_h5(self, input_particle,s=0):
        """
        Set the input particle distribution from an OpenPMD HDF5 object.

        If you are having a OpenPMD HDF5 file, please do

        input_particle = h5py.File(input_file, "r") first

        This Converts the HDF5 file to SDDS format using fbpic2sdds.

        Parameters
        ----------
        input_particle : h5py.object
            Read HDF5 object in OpenPMD format (e.g. from FBPIC output).
        s : starting point, int, optional
            Iteration index to read from the HDF5 file. Default is 0.

        Returns
        -------
        sdds.SDDS
            Converted SDDS object ready for Elegant input.
        """
        self.input = input_particle
        self.sdds_input = pyclara.Converters.fbpic2sdds(self.input, s)
        self.detect = True
        return self.sdds_input

    def set_input(self, input_particle, s=0):
        """
                Automatically detect input type and set the particle distribution.

                Calls set_input_sdds() for .sdds file paths and set_input_h5()
                for HDF5 file objects.

                Parameters
                ----------
                input_particle : str or h5py.File
                    Either a path to a .sdds file, or an open h5py.File object
                    in OpenPMD format.
                s : int, optional
                    Iteration index, only used for HDF5 input. Default is 0.

                Raises
                ------
                ValueError
                    If the input type is not a .sdds path or h5py.File object.
                """
        if isinstance(input_particle, str) and input_particle.endswith('.sdds'):

            self.set_input_sdds(input_particle)

        elif isinstance(input_particle, h5py.File):

            self.set_input_h5(input_particle, s)
        else:
            raise ValueError(f"Unrecognised input type: {type(input_particle)}. "
                             "Expected .sdds path or h5py.File object.")

    def write_elegant_input(self):
        """
                Write all Elegant input files to the run directory.

                Creates the run directory if it does not exist, then writes:
                - The lattice file ({name}.lte)
                - The run file ({name}.ele), with lattice, beamline, and input
                  beam updated to match the current simulation name and input
                - The input particle distribution ({name}_input.sdds)
                - Any wake files referenced in the lattice

                If the input was set from an HDF5 file (detect=True), Twiss
                parameters are also updated from the input distribution using
                fbpic2twiss. Otherwise the existing .ele file is used with
                string substitution.
        """
        os.makedirs(self.run_dir, exist_ok=True)
        filepath1 = os.path.join(self.run_dir, f'{self.name}.lte')
        filepath2 = os.path.join(self.run_dir, f'{self.name}.ele')
        filepath3 = os.path.join(self.run_dir, f'{self.name}_input.sdds' )
        with open(filepath1, 'w') as f:
            f.write(self.elegant_lte)

        if self.detect == True:
            f = pyclara.Converters.elegant_ele()
            f.load(self.ele_path)
            d = pyclara.Converters.fbpic2twiss(self.input)
            f.update_from_dict(d)
            f.run_setup["lattice"] = f'{self.name}.lte'
            f.run_setup["use_beamline"] = f'{self.name}'
            f.sdds_beam["input"] = f'{self.name}_input.sdds'
            f.write(filepath2)
        else:
            self.elegant_ele = self.elegant_ele.replace(f'{os.path.basename(self.sdds_input_path)}', f'{self.name}_input.sdds')
            self.elegant_ele = re.sub(r'lattice = \S+', f'lattice = {self.name}.lte' , self.elegant_ele)
            self.elegant_ele = re.sub(r'use_beamline = \S+', f'use_beamline = {self.name}', self.elegant_ele)
            with open(filepath2, 'w') as f:
                f.write(self.elegant_ele)

        self.sdds_input.mode = sdds.SDDS.SDDS_BINARY
        self.sdds_input.save(filepath3)
        for name, data in self.wakefiles.items():
            filepath = os.path.join(self.run_dir, f'{name}')
            with open(filepath, 'w') as f:
                f.write(data)

    def run (self):
        """
                Execute the Elegant simulation using MPI.

                Runs Pelegant with 6 MPI processes in the run directory.
                stdout and stderr are captured and printed.

                default Requires Pelegant to be installed at
                /elegant/elegant/bin/Linux-x86_64/Pelegant and mpirun at
                /lib64/mpich/bin/mpirun.

                But you always can set cmd as
                # e = Elegant_runner('/path/to/elegant_input', 'FEBE', '/path/to/output')
                # e.cmd = [
                #     'mpirun',       # or full path — find with: which mpirun
                #     '-np', '32',    # number of cores
                #     'Pelegant',     # or full path — find with: which Pelegant
                #     f'{e.name}.ele'
                # ]
                # e.run()
            """
        result = subprocess.run(
            self.cmd,
            cwd=self.run_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        print(result.stdout)

    def get_output(self):
        """
                Read the most recent Elegant output SDDS file and convert to HDF5.

                Finds the most recently modified .sdds file in the run directory
                and converts it to OpenPMD HDF5 format using sdds2fbpic.

                Returns
                -------
                h5py.File
                    Output particle distribution in OpenPMD HDF5 format.

                Raises
                ------
                ValueError
                    If no .sdds files are found in the run directory.
                """
        path = str(max(Path(self.run_dir).glob("*.sdds"), key=lambda f: f.stat().st_mtime))
        print(path)
        data = pyclara.Converters.sdds2fbpic(path)
        return data

    def track (self, particles):
        """
        Run a full Elegant tracking step in the S2E pipeline.

        Convenience method that calls set_input(), write_elegant_input(),
        run(), and get_output() in sequence. This is the method called
        by the S2E pipeline.

        Parameters
        ----------
        particles : str or h5py.File
            Input particle distribution. Either a path to a .sdds file
            or an open h5py.File object in OpenPMD format.

        Returns
        -------
        h5py.File
            Output particle distribution in OpenPMD HDF5 format, ready
            to be passed to the next tracker in the pipeline.

        Examples
        --------
        # output = e.track(h5py.File('input.h5', 'r'))
        """
        self.set_input(particles)
        self.write_elegant_input()
        self.run()
        output = self.get_output()
        return output









