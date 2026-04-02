import re as _re
from collections import defaultdict as _defaultdict
import sdds 
import numpy as np
import h5py
from scipy import constants

def elegant_lte_loader(filename):
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

        if etype.upper() == "LINE" :
            val = _re.findall(r'=\((.+)\)',params_str)
            elist = []
            for eval  in val[0].split(","):
                elist.append(eval.strip())
            params["LINE"] = elist

        params['TYPE'] = etype.upper()
        elements[name] = params

    return dict(elements)

def elegant_lte_writer(elements, filename):

    with open(filename, 'w') as f:

        # --- Part 1: write all element definitions ---
        for name, elem in elements.items():

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
            f.write(_elegant_longname_wrap(line) + '\n')

        # --- Part 2: write the LINE definition at the end ---
        line_key = next(k for k, v in elements.items() if v['TYPE'] == 'LINE')
        beamline  = elements[line_key]['LINE']

        chunks   = [beamline[i:i+6] for i in range(0, len(beamline), 6)]
        line_str = ', &\n'.join(', '.join(chunk) for chunk in chunks)
        f.write(f"\n{line_key}: LINE=({line_str})\n")


def _elegant_longname_wrap(line, width=100):

    if len(line) <= width:
        return line

    parts   = line.split(', ')
    result  = []
    current = ''

    for part in parts:
        if len(current) + len(part) > width:
            result.append(current.rstrip(', ') + ', &')
            current = part + ', '
        else:
            current += part + ', '

    result.append(current.rstrip(', '))
    return '\n'.join(result)

def elegant_lte_splitter(filename, split_element):

    # Step 1: load the full file
    elements = elegant_lte_loader(filename)

    # Step 2: find LINE key and the ordered element list
    line_key  = next(k for k, v in elements.items() if v['TYPE'] == 'LINE')
    full_line = elements[line_key]['LINE']

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

    elements2 = {'START': elements['START']} 

    for name, elem in elements.items():
        if elem['TYPE'] == 'LINE':
            continue
        if name.upper() in set1:
            elements1[name] = elem
        elif name.upper() in set2:
            elements2[name] = elem


    # Step 5: add LINE entries to each dict
    elements1[line_key + '_1'] = {'NAME': line_key + '_1', 'TYPE': 'LINE', 'LINE': line1}
    elements2[line_key + '_2'] = {'NAME': line_key + '_2', 'TYPE': 'LINE', 'LINE': line2}

    return elements1, elements2

class elegeant_ele:
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

    def read(self,filename):
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

    
    


def elegant_ele_writer(template_ele, output_ele, run_name,
                        lte_file, beamline, twiss=None, sdds_input=None):

    with open(template_ele, 'r') as f:
        content = f.read()

    def replace_param(text, key, value):
        return re.sub(rf'{key}\s*=\s*\S+', f'{key} = {value}', text)

    content = content.replace('%s', run_name)
    content = replace_param(content, 'lattice',      lte_file)
    content = replace_param(content, 'use_beamline', beamline)

    if twiss is not None:
        content = replace_param(content, 'p_central', twiss['p_central'])
        content = replace_param(content, 's_start',   twiss['s_start'])
        content = replace_param(content, 'Z0',        twiss['s_start'])
        content = replace_param(content, 'beta_x',    twiss['beta_x'])
        content = replace_param(content, 'alpha_x',   twiss['alpha_x'])
        content = replace_param(content, 'eta_x',     twiss['eta_x'])
        content = replace_param(content, 'etap_x',    twiss['etap_x'])
        content = replace_param(content, 'beta_y',    twiss['beta_y'])
        content = replace_param(content, 'alpha_y',   twiss['alpha_y'])
        content = replace_param(content, 'eta_y',     twiss['eta_y'])
        content = replace_param(content, 'etap_y',    twiss['etap_y'])

    if sdds_input is not None:
        content = replace_param(content, 'input', sdds_input)

    with open(output_ele, 'w') as f:
        f.write(content)


def fbpic2sdds(inputfile, outputfile, particles_group):

    x_all = []
    y_all= []
    z_all = []
    px_all = []
    py_all = []
    pz_all = []
    c_all = 0.0
    with h5py.File(inputfile, 'r') as f:

        iteration = list(f["data"].keys())[0]

        for particle in particles_group :
                x_all.append(f[f"/data/{iteration}/particles/{particle}/position/x"][:])
                y_all.append(f[f"/data/{iteration}/particles/{particle}/position/y"][:])
                z_all.append(f[f"/data/{iteration}/particles/{particle}/position/z"][:])
                px_all.append(f[f"/data/{iteration}/particles/{particle}/momentum/x"][:])
                py_all.append(f[f"/data/{iteration}/particles/{particle}/momentum/y"][:])
                pz_all.append(f[f"/data/{iteration}/particles/{particle}/momentum/z"][:])
                c_all += f[f"/data/{iteration}/particles/{particle}/charge"].attrs["value"]
    x  = np.concatenate(x_all)
    y  = np.concatenate(y_all)
    z  = np.concatenate(z_all)
    px = np.concatenate(px_all)
    py = np.concatenate(py_all)
    pz = np.concatenate(pz_all)
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
    sdds_obj.defineParameter("c", units='C', description='Charge of the beam')
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
    sdds_obj.setParameterValue("c", c_all , page=1)
    sdds_obj.setParameterValue("Particles", len(x) , page=1)
    for col_name, col_data in columns_page1.items():
        sdds_obj.setColumnValueList(col_name, col_data, page=1)

    sdds_obj.save(outputfile)

def sdds2fbpic(sddsfile) :
    f = sdds.load(sddsfile)
    xbeam = np.array(f.getColumnValueList("x"))
    ybeam = np.array(f.getColumnValueList("y"))
    xpbeam = np.array(f.getColumnValueList("xp"))
    ypbeam = np.array(f.getColumnValueList("yp"))
    dtbeam = np.array(f.getColumnValueList("dt"))
    zbeam = np.array(f.getColumnValueList("dt"))* constants.c
    pbeam = np.array(f.getColumnValueList("p"))
    pxbeam = xpbeam * pbeam
    pybeam = ypbeam * pbeam


    return {"x" : xbeam, "y" : ybeam, "z" : zbeam, "xp" : xpbeam, "yp" : ypbeam, "dt" : dtbeam, "p" : pbeam, "px" : pxbeam, "py" : pybeam}

def fbpic2twiss(inputfile , particles_group):
    x_all = []
    y_all= []
    z_all = []
    px_all = []
    py_all = []
    pz_all = []
    w_all = []
    c_all = 0.0
    with h5py.File(inputfile, 'r') as f:

        iteration = list(f["data"].keys())[0]

        for particle in particles_group :
                x_all.append(f[f"/data/{iteration}/particles/{particle}/position/x"][:])
                y_all.append(f[f"/data/{iteration}/particles/{particle}/position/y"][:])
                z_all.append(f[f"/data/{iteration}/particles/{particle}/position/z"][:])
                px_all.append(f[f"/data/{iteration}/particles/{particle}/momentum/x"][:])
                py_all.append(f[f"/data/{iteration}/particles/{particle}/momentum/y"][:])
                pz_all.append(f[f"/data/{iteration}/particles/{particle}/momentum/z"][:])
                w_all.append(f[f"/data/{iteration}/particles/{particle}/weighting"][:])
                c_all += f[f"/data/{iteration}/particles/{particle}/charge"].attrs["value"]
    x  = np.concatenate(x_all)
    y  = np.concatenate(y_all)
    z  = np.concatenate(z_all)
    px = np.concatenate(px_all)
    py = np.concatenate(py_all)
    pz = np.concatenate(pz_all)
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