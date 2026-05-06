import re as _re
from collections import defaultdict as _defaultdict
import sdds 
import numpy as np
import h5py
from scipy import constants


def elegant_lte_loader(filename):
    """
    Load an Elegant lattice file (.lte) into a dictionary.

    Parameters
    ----------
    filename : str
        Path to the Elegant .lte lattice file.

    Returns
    -------
    dict
        Dictionary of {element_name: params_dict} where each
        params_dict contains keys such as 'NAME', 'TYPE', and
        any element-specific parameters (e.g. 'L', 'K1', 'WAKEFILE').
        LINE elements additionally contain a 'LINE' key with an
        ordered list of element names.
    """
    with open(filename) as f:
        content = f.read()

    # Remove comments (! and %)
    content = _re.sub(r'[!%].*', '', content)

    # Join continuation lines
    content = _re.sub(r'&\s*\n', ' ', content)

    # Find all name&type blocks: NAME: TYPE, key=val, key=val ...
    pattern = _re.compile(
        r'([\w-]+)\s*:\s*(\w+)\s*,?\s*([^;]*)',
        _re.IGNORECASE
    )

    elements = _defaultdict(list)

    for match in pattern.finditer(content):
        name, etype, params_str = match.groups()
        params = {"NAME": name.upper()}

        for kv in _re.findall(r'(\w+)\s*=\s*([^,\n]+)', params_str):
            key, val = kv
            val = val.strip().strip('"')
            try:
                params[key.upper()] = float(val)
            except ValueError:
                params[key.upper()] = val

        if etype.upper() == "LINE":
            val = _re.findall(r'=\((.+)\)', params_str)
            elist = []
            for eval in val[0].split(","):
                elist.append(eval.strip())
            params["LINE"] = elist

        params['TYPE'] = etype.upper()
        elements[name] = params


    return dict(elements)

class elegant_lte:
    """
    Open an Elegant lattice file (.lte) as a Python object.

    Provides methods to load, split, and write Elegant lattice files.
    Internally stores all elements as a dictionary.

    Parameters
    ----------
    filename : str, optional
        Path to an Elegant .lte file to load on initialisation.
    elements : dict, optional
        Pre-built elements dictionary (e.g. from splitter()).
        If provided, the object is initialised directly from this
        dict without reading a file.

    Attributes
    ----------
    filename : str
        Path to the loaded lattice file.
    elements_dict : dict
        Dictionary of {element_name: params_dict} for all elements
        in the lattice.

    Notes
    -----
    A typical use case for the S2E pipeline is to split the full
    CLARA-FEBE lattice into two sections — before and after the
    plasma stage — so that each section can be tracked separately
    by TrackerElegant, with a plasma simulation in between.

    Examples
    --------
    Split the FEBE lattice at the plasma entrance element and write
    two sub-lattice files for use in the S2E pipeline:

    import pyclara
    inputdir = '/path/to/elegant_input'
    f = pyclara.Converters.elegant_lte(filename=f'{inputdir}/FEBE.lte')
    f.load()
    # Split at the plasma entrance element
    f1, f2 = f.splitter('CLA-FEC1-SIM-FOCUS-01')
    # Write the two sub-lattices to disk
    f1.writer(f'{inputdir}/FEBE1.lte')   # before plasma
    f2.writer(f'{inputdir}/FEBE2.lte')   # after plasma
    """
    def __init__(self, filename = None, elements = None):
        if filename is not None:
            self.filename = filename
            self.load()
        if elements is not None:
            self.elements_dict = elements

    def load(self):
        """
        Load an Elegant lattice file (.lte) into a dictionary.

        Parameters
        ----------
        filename : str
            Path to the Elegant .lte lattice file.

        Returns
        -------
        dict
            Dictionary of {element_name: params_dict} where each
            params_dict contains keys such as 'NAME', 'TYPE', and
            any element-specific parameters (e.g. 'L', 'K1', 'WAKEFILE').
            LINE elements additionally contain a 'LINE' key with an
            ordered list of element names.
    """
        with open(self.filename) as f:
            content = f.read()

        # Remove comments (! and %)
        content = _re.sub(r'[!%].*', '', content)

        # Join continuation lines
        content = _re.sub(r'&\s*\n', ' ', content)

        # Find all name&type blocks: NAME: TYPE, key=val, key=val ...
        pattern = _re.compile(
            r'([\w-]+)\s*:\s*(\w+)\s*,?\s*([^;]*)',
            _re.IGNORECASE
        )

        elements = _defaultdict(list)

        for match in pattern.finditer(content):
            name, etype, params_str = match.groups()
            params = {"NAME": name.upper()}

            for kv in _re.findall(r'(\w+)\s*=\s*([^,\n]+)', params_str):
                key, val = kv
                val = val.strip().strip('"')
                try:
                    params[key.upper()] = float(val)
                except ValueError:
                    params[key.upper()] = val

            if etype.upper() == "LINE" :
                val = _re.findall(r'=\((.+)\)',params_str)
                elist = []
                for eval  in val[0].split(","):
                    elist.append(eval.strip())
                params["LINE"] = elist

            params['TYPE'] = etype.upper()
            elements[name] = params
            self.elements_dict =  dict(elements)

        return self.elements_dict

    def splitter(self, split_element):
        """
        Split the beamline at a named element into two sub-lattices.

        Divides the full beamline into two elegant_lte objects:
        the first contains all elements up to and including
        split_element, the second contains all elements after it.

        Parameters
        ----------
        split_element : str
            Name of the element at which to split the beamline.
            Must match exactly (case-sensitive).

        Returns
        -------
        lte1 : elegant_lte
            Sub-lattice from the start up to and including
            split_element.
        lte2 : elegant_lte
            Sub-lattice from after split_element to the end,
            with a 'START' marker prepended.

        """

        line_key  = next(k for k, v in self.elements_dict.items() if v['TYPE'] == 'LINE')
        full_line = self.elements_dict[line_key]['LINE']

        # Step 3: find split index — exact match
        if split_element not in full_line:

        # check if a case-insensitive match exists
            close_match = next(
                (e for e in full_line if e.upper() == split_element.upper()),
                None
            )

            if close_match:
                raise ValueError(
                    f"'{split_element}' not found. Did you mean '{close_match}'?"
                )
            else:
                raise ValueError(
                f"'{split_element}' not found in beamline."
                )

        split_idx = full_line.index(split_element)


        line1 = full_line[:split_idx+1]    # before split element
        line2 = ['START'] + full_line[split_idx+1:]    # from split element onwards

        # Step 4: build two element dicts
        set1 = set(e.upper() for e in line1)
        set2 = set(e.upper() for e in line2)

        elements1 = {}
        elements2 = {}

        elements2 = {'START': self.elements_dict['START']}

        for name, elem in self.elements_dict.items():
            if elem['TYPE'] == 'LINE':
                continue
            if name.upper() in set1:
                elements1[name] = elem
            elif name.upper() in set2:
                elements2[name] = elem


        # Step 5: add LINE entries to each dict
        elements1[line_key + '1'] = {'NAME': line_key + '_1', 'TYPE': 'LINE', 'LINE': line1}
        elements2[line_key + '2'] = {'NAME': line_key + '_2', 'TYPE': 'LINE', 'LINE': line2}

        return elegant_lte(elements=elements1), elegant_lte(elements=elements2)

    def writer(self, outputfile):
        """
        Write the lattice to an Elegant .lte file.

        Writes all element definitions followed by the LINE
        definition. Long lines are automatically wrapped with
        Elegant continuation syntax (', &').

        Parameters
        ----------
        outputfile : str
            Path to the output .lte file to write.
        """
        with open(outputfile, 'w') as f:

            # --- Part 1: write all element definitions ---
            for name, elem in self.elements_dict.items():

                if elem['TYPE'] == 'LINE':
                    continue

                params = {k: v for k, v in elem.items() if k not in ('NAME', 'TYPE')}

                param_parts = []
                for k, v in params.items():
                    if isinstance(v, float):
                        param_parts.append(f'{k} = {v:g}')
                    else:
                        param_parts.append(f'{k} = "{v}"')

                line = f"{name}: {elem['TYPE']}, " + ', '.join(param_parts) + ';'
                f.write(self._longname_wrap(line) + '\n')

            # --- Part 2: write the LINE definition at the end ---
            line_key = next(k for k, v in self.elements_dict.items() if v['TYPE'] == 'LINE')
            beamline = self.elements_dict[line_key]['LINE']

            chunks = [beamline[i:i + 6] for i in range(0, len(beamline), 6)]
            line_str = ', &\n'.join(', '.join(chunk) for chunk in chunks)
            f.write(f"\n{line_key}: LINE=({line_str})\n")

    def _longname_wrap(self, line, width=100):

        if len(line) <= width:
            return line

        parts = line.split(', ')
        result = []
        current = ''

        for part in parts:
            if len(current) + len(part) > width:
                result.append(current.rstrip(', ') + ', &')
                current = part + ', '
            else:
                current += part + ', '

        result.append(current.rstrip(', '))
        return '\n'.join(result)

class elegant_ele:
    """
        Represents an Elegant run file (.ele) as a Python object.

        Provides methods to load, update, and write Elegant .ele files.
        Each section of the .ele file (run_setup, sdds_beam, etc.) is
        stored as a dictionary attribute.

        Attributes
        ----------
        global_settings : dict
        Controls global Elegant settings such as inhibit_fsync.
        Corresponds to the &global_settings block in the .ele file.

        run_setup : dict
            Core simulation setup. Must include:
            - 'lattice'       : path to the .lte lattice file
            - 'use_beamline'  : name of the beamline to track through
            - 'p_central'     : reference momentum βγ [m_e*c]
            - 's_start'       : starting s position along the beamline [m]
            Updated automatically by update_from_dict() after a plasma stage.
            Corresponds to the &run_setup block.

        run_control : dict
            Controls the number of simulation passes.
            Key entry: 'n_passes' (integer, usually 1 for single-pass tracking).
            Corresponds to the &run_control block.

        twiss_output : dict
            Initial Twiss parameters at the start of the beamline.
            Key entries: 'beta_x', 'beta_y', 'alpha_x', 'alpha_y',
            'eta_x', 'eta_y', 'etap_x', 'etap_y'.
            Updated automatically by update_from_dict() when chaining
            after a FBPIC plasma stage, using values computed by fbpic2twiss().
            Corresponds to the &twiss_output block.

        floor_coordinates : dict
            Physical floor coordinates of the beamline start point.
            Key entry: 'Z0' (longitudinal position [m]), updated automatically
            by update_from_dict() to match s_start after a plasma stage.
            Corresponds to the &floor_coordinates block.

        matrix_output : dict
            Controls transfer matrix output along the beamline.
            Corresponds to the &matrix_output block.

        sdds_beam : dict
            Specifies the input beam SDDS file.
            Key entry: 'input' (filename of the .sdds input beam file).
            Set automatically by write_elegant_input() in TrackerElegant
            to '{name}_input.sdds'.
            Corresponds to the &sdds_beam block.

        track : dict
            Executes the tracking. Typically left empty — the presence
            of the &track block is sufficient to trigger tracking in Elegant.
            Corresponds to the &track block.


    """
    def __init__(self):
        self.global_settings = {}
        self.run_setup = {}
        self.run_control = {}
        self.twiss_output = {}
        self.floor_coordinates = {}
        self.matrix_output = {}
        self.sdds_beam = {}
        self.track = {}

        #self.global_settings["inhibit_fsync"] = 0

    def write(self,filename):
        """
        Write the current state to an Elegant .ele file.

        Parameters
        ----------
        filename : str
            Path to the output .ele file to write.

        Examples
        --------
        ele.write('/path/to/run/FEBE.ele')
        """
        f = open(filename, "w")
        self._writedict(f,"global_settings",self.global_settings)
        self._writedict(f,"run_setup",self.run_setup)
        self._writedict(f,"run_control",self.run_control)
        self._writedict(f,"twiss_output",self.twiss_output)
        self._writedict(f,"floor_coordinates",self.floor_coordinates)
        self._writedict(f,"matrix_output",self.matrix_output)
        self._writedict(f,"sdds_beam",self.sdds_beam)
        self._writedict(f,"track",self.track)
        f.close()

    def _writedict(self, f, name, d):
        f.write(f'&{name}\n')
        for k, v in d.items() :
            f.write(f'\t{k} = {v}\n')
        f.write(f'&end\n')

    def load(self,filename):
        """
        Load and parse an existing Elegant .ele file.

        Reads all recognised &block ... &end sections and populates
        the corresponding dictionary attributes. Unrecognised blocks
        are silently ignored.

        Parameters
        ----------
        filename : str
            Path to the Elegant .ele file to load.
        """
        f = open(filename, "r")
        global_settings = False
        run_control     = False
        run_setup       = False
        twiss_output    = False
        floor_coordinates = False
        matrix_output   = False
        sdds_beam       = False
        track           = False
        
        for l in f :
            if l.strip() == "&run_control":
                run_control = True
                continue
            elif l.strip() == "&global_settings":
                global_settings = True
                continue
            elif l.strip() == "&run_setup":
                run_setup = True
                continue
            elif l.strip() == "&twiss_output":
                twiss_output = True
                continue
            elif l.strip() == "&floor_coordinates":
                floor_coordinates = True
                continue
            elif l.strip() == "&matrix_output":
                matrix_output = True
                continue
            elif l.strip() == "&sdds_beam":
                sdds_beam = True
                continue
            elif l.strip() == "&track":
                track = True
                continue
                
            if l.strip() == "&end":
                if run_control: 
                    run_control = False
                elif global_settings:
                    global_settings = False
                elif run_setup:
                    run_setup = False   
                elif twiss_output:
                    twiss_output = False
                elif floor_coordinates:
                    floor_coordinates = False
                elif matrix_output:
                    matrix_output = False
                elif sdds_beam:
                    sdds_beam = False
                elif track:
                    track = False
                continue
            
            if run_control:
                split_line = l.split()
                self.run_control[split_line[0]] = split_line[2]
            if global_settings:
                split_line = l.split()
                self.global_settings[split_line[0]] = split_line[2]
            if run_setup:
                split_line = l.split()
                self.run_setup[split_line[0]] = split_line[2]
            if twiss_output:
                split_line = l.split()
                self.twiss_output[split_line[0]] = split_line[2]
            if floor_coordinates:
                split_line = l.split()
                self.floor_coordinates[split_line[0]] = split_line[2]
            if matrix_output:
                split_line = l.split()
                self.matrix_output[split_line[0]] = split_line[2]
            if sdds_beam:
                split_line = l.split()
                self.sdds_beam[split_line[0]] = split_line[2]
            if track:
                split_line = l.split()
                self.track[split_line[0]] = split_line[2]

        f.close()

    def update_from_dict(self, d):
        """
                Update run parameters from a dictionary.

                Typically called with the output of fbpic2twiss() to update
                the reference momentum, s-start position, and Twiss parameters
                after a FBPIC plasma stage, ready for the next Elegant section.
        +------------+---------------------+---------------------------+
        | Key        | Elegant block       | Description               |
        +============+=====================+===========================+
        | p_central  | run_setup           | Reference momentum βγ     |
        | s_start    | run_setup           | Starting s position [m]   |
        | s_start    | floor_coordinates   | Also sets Z0 = s_start    |
        | beta_x     | twiss_output        | Horizontal beta [m]       |
        | beta_y     | twiss_output        | Vertical beta [m]         |
        | alpha_x    | twiss_output        | Horizontal alpha          |
        | alpha_y    | twiss_output        | Vertical alpha            |
        | eta_x      | twiss_output        | Horizontal dispersion [m] |
        | eta_y      | twiss_output        | Vertical dispersion [m]   |
        | etap_x     | twiss_output        | Horizontal disp. prime    |
        | etap_y     | twiss_output        | Vertical disp. prime      |
        | input      | sdds_beam           | Input SDDS filename       |
         +------------+---------------------+---------------------------+

        Parameters
        ----------
        d : dict
            Dictionary of parameter updates, typically from fbpic2twiss().
            Unrecognised keys are silently ignored.

        """
        for k, v in d.items():
            if k == "p_central":
                self.run_setup["p_central"] = v
            if k == "s_start":
                self.run_setup["s_start"] = v
                self.floor_coordinates['Z0'] = v
            if k == 'beta_x':
                self.twiss_output["beta_x"] = v
            if k == 'beta_y':
                self.twiss_output["beta_y"] = v
            if k == 'alpha_x':
                self.twiss_output["alpha_x"] = v
            if k == 'alpha_y':
                self.twiss_output["alpha_y"] = v
            if k == 'eta_x':
                self.twiss_output["eta_x"] = v
            if k == 'eta_y':
                self.twiss_output["eta_y"] = v
            if k == 'etap_x':
                self.twiss_output["etap_x"] = v
            if k == 'etap_y':
                self.twiss_output["etap_y"] = v
            if k == 'input':
                self.sdds_beam['input'] = v

def fbpic2sdds(input, z_offset = 0,  particles_group = ["beam"]):
    """
        Convert an FBPIC OpenPMD HDF5 output to an Elegant SDDS object.

        Reads particle positions and momenta from an OpenPMD HDF5 file,
        converts from FBPIC coordinates (x, y, z, ux, uy, uz) to Elegant
        coordinates (x, xp, y, yp, t, p, dt), and returns an SDDS object
        ready to be written and used as Elegant beam input.

        Coordinate conventions
        ----------------------
        - FBPIC momenta (ux, uy, uz) are in SI units [kg·m/s], normalised
          to m_e*c to give βγ.
        - Elegant uses xp = px/p, yp = py/p (divergence angles) and
          p = |βγ| (total momentum in units of m_e*c).
        - Time coordinate t = z/c, relative time dt = t - mean(t).

        Parameters
        ----------
        input : h5py.File
            Open HDF5 file object in OpenPMD format from FBPIC output.
        z_offset : float, optional
            Longitudinal offset to add to z positions [m]. Default is 0.
        particles_group : list of str, optional
            List of particle species names to include. 'electrons' is
            excluded (background plasma). Default is ['beam'].

        Returns
        -------
        sdds.SDDS
            SDDS object containing columns x, xp, y, yp, t, p, dt,
            particleID and parameters pCentral, charge, Particles.

        Examples
        --------
        import h5py
        f = h5py.File('fbpic_output.h5', 'r')
        sdds_obj = fbpic2sdds(f, z_offset=0.1)
        sdds_obj.save('beam_input.sdds')
        """
    x_all = []
    y_all= []
    z_all = []
    px_all = []
    py_all = []
    pz_all = []
    c_all = 0.0

    iteration = list(input["data"].keys())[0]

    for particle in particles_group:
        if particle != "electrons":
            x_all.append(input[f"/data/{iteration}/particles/{particle}/position/x"][:])
            y_all.append(input[f"/data/{iteration}/particles/{particle}/position/y"][:])
            z_all.append(input[f"/data/{iteration}/particles/{particle}/position/z"][:])
            px_all.append(input[f"/data/{iteration}/particles/{particle}/momentum/x"][:])
            py_all.append(input[f"/data/{iteration}/particles/{particle}/momentum/y"][:])
            pz_all.append(input[f"/data/{iteration}/particles/{particle}/momentum/z"][:])
            c_all += input[f"/data/{iteration}/particles/{particle}/charge"].attrs["value"]
    x  = np.concatenate(x_all)
    y  = np.concatenate(y_all)
    z  = np.concatenate(z_all) + z_offset
    px = np.concatenate(px_all)/constants.m_e/constants.c
    py = np.concatenate(py_all)/constants.m_e/constants.c
    pz = np.concatenate(pz_all)/constants.m_e/constants.c
    n  = np.arange(1, len(x)+1)

    t = z / constants.c
    dt = t - np.mean(t)
    p = np.sqrt(px**2+py**2+pz**2)
    xp = px/p
    yp = py/p


    sdds_obj = sdds.SDDS()
    
    sdds_obj.setDescription(
        "SDDS files converted from hdf5 files from FBPIC",
        "x xp y yp t p dt" 
    )
    sdds_obj.defineParameter("pCentral", symbol="p$bcen$n", units='m$be$nc', description="Reference beta*gamma",type=sdds.SDDS_DOUBLE)
    sdds_obj.defineParameter("charge", units='C', description="Bunch charge before sampling")
    sdds_obj.defineParameter("Particles", description='Number of particles',type=sdds.SDDS_LONG)
    sdds_obj.defineColumn("x",  units="m")
    sdds_obj.defineColumn("xp")
    sdds_obj.defineColumn("y",  units="m")
    sdds_obj.defineColumn("yp")
    sdds_obj.defineColumn("t",  units="s")
    sdds_obj.defineColumn("p",  units="m$be$nc")
    sdds_obj.defineColumn("dt", units="s")
    sdds_obj.defineColumn("particleID",type=sdds.SDDS_LONG)
    columns_page1 = {
        "x"  : x.tolist(),
        "xp" : xp.tolist(),
        "y"  : y.tolist(),
        "yp" : yp.tolist(),
        "t"  : t.tolist(),
        "p"  : p.tolist(),
        "dt" : dt.tolist(),
        "particleID" : n.tolist()
    }
    sdds_obj.setParameterValue("charge", c_all , page=1)
    sdds_obj.setParameterValue("Particles", len(x) , page=1)
    sdds_obj.setParameterValue("pCentral", np.average(p), page=1)
    for col_name, col_data in columns_page1.items():
        sdds_obj.setColumnValueList(col_name, col_data, page=1)

    sdds_obj.mode = sdds.SDDS.SDDS_BINARY
    return sdds_obj

def sdds2fbpic(sddsfile) :
    """
    Convert an Elegant SDDS output file to a particle dictionary.

    Reads particle coordinates from an Elegant SDDS output file and
    returns them as a dictionary of numpy arrays in FBPIC-compatible
    units.

    Coordinate conventions
    ----------------------
    - z is computed as t * c (longitudinal position from time coordinate).
    - px = xp * p, py = yp * p (momenta reconstructed from divergence
      angles and total momentum βγ).

    Parameters
    ----------
    sddsfile : str
        Path to the Elegant output SDDS file.

    Returns
    -------
    dict
        Dictionary with keys:
        'x', 'y'   — transverse positions [m]
        'z'        — longitudinal position [m]
        'xp', 'yp' — transverse divergence angles [rad]
        'dt'       — relative time offset [s]
        'p'        — total momentum βγ [dimensionless]
        'px', 'py' — transverse momenta βγ [dimensionless]

    Examples
    --------
    d = sdds2fbpic('FEBE.sdds')
    x = d['x']
    pz = d['p']
    """
    f = sdds.load(sddsfile)
    xbeam = np.array(f.getColumnValueList("x"))
    ybeam = np.array(f.getColumnValueList("y"))
    xpbeam = np.array(f.getColumnValueList("xp"))
    ypbeam = np.array(f.getColumnValueList("yp"))
    dtbeam = np.array(f.getColumnValueList("dt"))
    zbeam = np.array(f.getColumnValueList("t"))* constants.c #3change has been made in here dt
    pbeam = np.array(f.getColumnValueList("p"))
    pxbeam = xpbeam * pbeam
    pybeam = ypbeam * pbeam


    return {"x" : xbeam, "y" : ybeam, "z" : zbeam, "xp" : xpbeam, "yp" : ypbeam, "dt" : dtbeam, "p" : pbeam, "px" : pxbeam, "py" : pybeam}

def fbpic2twiss(input, particles_group = ['beam']):
    """
        Compute Twiss parameters and dispersion from an FBPIC HDF5 output.

        Reads particle data from an OpenPMD HDF5 file and computes the
        weighted second-moment Twiss parameters (beta, alpha), dispersion
        (eta, etap), and reference momentum for use as initial conditions
        in the subsequent Elegant tracking section.

        Parameters
        ----------
        input : h5py.File
            Open HDF5 file object in OpenPMD format from FBPIC output.
        particles_group : list of str, optional
            List of particle species to include in the calculation.
            Default is ['beam'].

        Returns
        -------
        dict
            Dictionary with keys:
            'p_central' — weighted mean longitudinal momentum βγ
            's_start'   — weighted mean z position [m]
            'beta_x'    — horizontal beta function [m]
            'alpha_x'   — horizontal alpha (=-0.5 * dbeta/ds)
            'eta_x'     — horizontal dispersion [m]
            'etap_x'    — horizontal dispersion prime
            'beta_y'    — vertical beta function [m]
            'alpha_y'   — vertical alpha
            'eta_y'     — vertical dispersion [m]
            'etap_y'    — vertical dispersion prime

        Notes
        -----
        Twiss parameters are computed from the beam sigma matrix:
            eps_x  = sqrt(<x²><x'²> - <xx'>²)
            beta_x = <x²> / eps_x
            alpha_x = -<xx'> / eps_x

        Dispersion is computed as:
            eta_x = <x·delta> / <delta²>

        where delta = (pz - p_central) / p_central.

        Examples
        --------
        d = fbpic2twiss(h5py.File('fbpic_output.h5', 'r'))
        print(d['beta_x'], d['alpha_x'])
        """
    x_all = []
    y_all= []
    z_all = []
    px_all = []
    py_all = []
    pz_all = []
    w_all = []
    c_all = 0.0

    iteration = list(input["data"].keys())[0]

    for particle in particles_group :
            x_all.append(input[f"/data/{iteration}/particles/{particle}/position/x"][:])
            y_all.append(input[f"/data/{iteration}/particles/{particle}/position/y"][:])
            z_all.append(input[f"/data/{iteration}/particles/{particle}/position/z"][:])
            px_all.append(input[f"/data/{iteration}/particles/{particle}/momentum/x"][:])
            py_all.append(input[f"/data/{iteration}/particles/{particle}/momentum/y"][:])
            pz_all.append(input[f"/data/{iteration}/particles/{particle}/momentum/z"][:])
            w_all.append(input[f"/data/{iteration}/particles/{particle}/weighting"][:])
            c_all += input[f"/data/{iteration}/particles/{particle}/charge"].attrs["value"]
    x  = np.concatenate(x_all)
    y  = np.concatenate(y_all)
    z  = np.concatenate(z_all)
    px = np.concatenate(px_all)/constants.m_e/constants.c
    py = np.concatenate(py_all)/constants.m_e/constants.c
    pz = np.concatenate(pz_all)/constants.m_e/constants.c
    w = np.concatenate(w_all)

    p_central = np.average(pz, weights=w)
    z_mean = np.average(z, weights=w)
    xp = px / pz
    yp = py / pz
    delta = (pz - p_central) / p_central

    # second moments
    x2  = np.average(x**2,  weights=w)
    xp2 = np.average(xp**2, weights=w)
    xxp = np.average(x*xp,  weights=w)
    y2  = np.average(y**2,  weights=w)
    yp2 = np.average(yp**2, weights=w)
    yyp = np.average(y*yp,  weights=w)

    # emittance
    eps_x = np.sqrt(x2*xp2 - xxp**2)
    eps_y = np.sqrt(y2*yp2 - yyp**2)

    # Twiss
    beta_x  =  x2  / eps_x
    alpha_x = -xxp / eps_x
    beta_y  =  y2  / eps_y
    alpha_y = -yyp / eps_y

    # dispersion
    d2     = np.average(delta**2, weights=w)
    eta_x  = np.average(x  * delta, weights=w) / d2
    etap_x = np.average(xp * delta, weights=w) / d2
    eta_y  = np.average(y  * delta, weights=w) / d2
    etap_y = np.average(yp * delta, weights=w) / d2

    return {
        'p_central': p_central,
        's_start' : z_mean,
        'beta_x':  beta_x,  'alpha_x': alpha_x,
        'eta_x':   eta_x,   'etap_x':  etap_x,
        'beta_y':  beta_y,  'alpha_y': alpha_y,
        'eta_y':   eta_y,   'etap_y':  etap_y,
    }

def sdds2beam_dict(sddsfile) :
    """
        Convert an Elegant SDDS file to a centred beam dictionary.

        Reads particle coordinates from an SDDS file, removes the mean
        from each coordinate (centres the distribution), and returns
        normalised phase space variables suitable for emittance analysis
        or FBPIC injection.

        Parameters
        ----------
        sddsfile : str
            Path to the Elegant SDDS file.

        Returns
        -------
        dict
            Dictionary with keys:
            'x'     — centred horizontal position [m]
            'xp'    — centred horizontal divergence [rad]
            'y'     — centred vertical position [m]
            'yp'    — centred vertical divergence [rad]
            'zeta'  — longitudinal position zeta = t*c [m]
            'delta' — relative momentum deviation (p - p0) / p0

        Notes
        -----
        Unlike sdds2fbpic(), this function centres all coordinates
        and returns delta (relative energy deviation) rather than
        absolute momentum. Useful for Twiss analysis and emittance
        calculations.

        Examples
        --------
        d = sdds2beam_dict('FEBE_output.sdds')
        print(d['delta'].std())   # RMS energy spread
        """
    f = sdds.load(sddsfile)
    x = np.array(f.getColumnValueList("x"))
    y = np.array(f.getColumnValueList("y"))
    xp = np.array(f.getColumnValueList("xp"))
    yp = np.array(f.getColumnValueList("yp"))
    t = np.array(f.getColumnValueList("t"))
    p = np.array(f.getColumnValueList("p"))

    x -= x.mean()
    y -= y.mean()
    xp -= xp.mean()
    yp -= yp.mean()
    t -= t.mean()
    p0 = p.mean()
    p -= p0

    zeta = t*constants.c
    delta = p/p0

    return {"x" : x,
            "xp" : xp,
            "y" : y,
            "yp" : yp,
            "zeta" : zeta,
            "delta" : delta}