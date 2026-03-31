import numpy as _np
from scipy import constants as _constants

m_e_eV = _constants.m_e*_constants.c**2/_constants.electron_volt

import sdds as _sdds
from ._elegant import elegant_lte_loader as _elegant_lte_loader

try :
    import xsuite as _xsuite
    import xobjects as _xobjects
    import xtrack as _xtrack
except ImportError:
    print("xsuite not found, cannot convert yaml to xsuite line")

def elegant2xsuite(elegant_file,
                   line_name = "FEBE",
                   start_element="CLA-FEA-SIM-DIP-04-END",
                   end_element="CLA-FED-SIM-DUMP-01-START",
                   elegant_twi = None,
                   elegant_ps = None) :

    # load elegant lattice file
    if isinstance(elegant_file, str) :
        lte = _elegant_lte_loader(elegant_file)

    # set default start/end
    if start_element == "" or start_element is None :
        lte_keys_list = list(lte.keys())
        print(lte_keys_list[0])
        start_element = lte[lte_keys_list[0]]['NAME']

    if end_element == "" or end_element is None :
        lte_keys_list = list(lte.keys())
        print(lte_keys_list[-2])
        end_element = lte[lte_keys_list[-1]]['NAME']

    # load elegant twiss file if provided
    if elegant_twi is not None :
        if isinstance(elegant_twi, str) :
            elegant_twi = _sdds.load(elegant_twi)

    # load elegant phase space file if provided
    if elegant_ps is not None :
        if isinstance(elegant_ps, str) :
            elegant_ps = _sdds.load(elegant_ps)

    # create xsuite environment
    env = _xtrack.Environment()

    # loop over elements
    for k in lte :
        ee = lte[k] # elegant element

        # skip over LINE
        if ee['TYPE'].upper() == "LINE":
            continue

        # Need to replace any - as xsuite will try to perform a deferred calculation
        name = ee['NAME'].replace("-","_")

        xe = None
        if ee['TYPE'] == 'CHARGE' :
            env.new(name, _xtrack.Marker)
        elif ee['TYPE'] == 'DRIFT':
            if ee['L'] == 0 :
                env.new(name,_xtrack.Marker)
            else :
                env[name+".L"] = float(ee['L'])
                env.new(name, _xtrack.Drift, length = name + ".L")
        elif ee['TYPE'] == 'CSRDRIFT':
            env[name+".L"] = float(ee['L'])
            env.new(name, _xtrack.Drift, length=name + ".L")
        elif ee['TYPE'] == 'LSCDRIFT':
            env[name+".L"] = float(ee['L'])
            env.new(name, _xtrack.Drift, length=name + ".L")
        elif ee['TYPE'] == "CSRCSBEND" or ee['TYPE'] == "CSBEND" :
            env[name+".L"] = float(ee['L'])
            env[name+".ANGLE"] = float(ee['ANGLE'])
            env[name+".E1"] = float(ee['E1'])
            env[name+".E2"] = float(ee['E2'])
            env.new(name, _xtrack.Bend,
                    length = name+".L",
                    angle = name+".ANGLE",
                    edge_entry_angle = name+".E1",
                    edge_exit_angle = name+".E2")
        elif ee['TYPE'] == 'KQUAD':
            env[name+".L"] = float(ee['L'])
            env[name+".K1"] = float(ee['K1'])
            env.new(name, _xtrack.Quadrupole, length=name+".L", k1=name+".K1")
        elif ee['TYPE'] == "KSEXT":
            env[name+".L"] = float(ee['L'])
            env[name+".K2"] = float(ee['K2'])
            env.new(name, _xtrack.Sextupole, length=name+".L", k2=name+".K2")
        elif ee['TYPE'] == "KICKER":
            env[name + ".L"] = float(ee['L'])
            env.new(name, _xtrack.Drift, length=name + ".L")
        elif ee['TYPE'] == "ECOL":
            env[name + ".L"] = float(ee['L'])
            env.new(name, _xtrack.Drift, length=name + ".L")
        elif ee['TYPE'] == "MAXAMP" :
            env.new(name, _xtrack.Marker)
        elif ee['TYPE'] == 'WATCH':
            env.new(name,_xtrack.Marker)
        elif ee['TYPE'] == 'MONI' :
            env[name + ".L"] = float(ee['L'])
            env.new(name, _xtrack.Drift, length=name+".L")
        elif ee['TYPE'] == "RFCW" :
            env[name+".L"] = float(ee['L'])
            env[name+".VOLT"] = float(ee['VOLT'])
            env[name+".FREQ"] = float(ee["FREQ"])
            env[name+".PHAS"] = float(ee["PHASE"])
            env.new(name, _xtrack.Cavity,
                    length=name+".L",
                    voltage=name+".VOLT",
                    frequency=name+".FREQ",
                    lag=name+".PHAS")
        elif ee['TYPE'] == "RFDF" :
            env[name+".L"] = float(ee['L'])
            env.new(name, _xtrack.Drift, length=name+".L")
        else :
            print("element type ", ee['TYPE'], " not recognised, skipping")

    # need to be uppercase
    line_component_names = [e.upper() for e in lte[line_name]['LINE']]

    if start_element is None or start_element == "":
        start_element = line_component_names[0]
    if end_element is None or end_element == "":
        end_element = line_component_names[-1]

    # select components within range
    subline_component_names = []

    adding = False
    istart = -1

    for i, e in enumerate(line_component_names) :
        if e == start_element :
            adding = True
            istart = i

        if adding :
            subline_component_names.append(e.replace("-","_"))

        if e == end_element :
            adding = False

    # create line
    xtrack_line = env.new_line(name = line_name,
                               components = subline_component_names)
    xtrack_twiss0 = None
    dict_twiss0 = None
    xtrack_twiss = None
    xtrack_particles = None

    # add twiss to environment if provided
    if elegant_twi is not None and elegant_ps is None:
        s = elegant_twi.getColumnValueList('s')[istart]
        p0 = elegant_twi.getColumnValueList('pCentral0')[istart]
        betax = elegant_twi.getColumnValueList('betax')[istart]
        alphax = elegant_twi.getColumnValueList('alphax')[istart]
        etax = elegant_twi.getColumnValueList('etax')[istart]
        etaxp = elegant_twi.getColumnValueList('etaxp')[istart]

        betay = elegant_twi.getColumnValueList('betay')[istart]
        alphay = elegant_twi.getColumnValueList('alphay')[istart]
        etay = elegant_twi.getColumnValueList('etay')[istart]
        etayp = elegant_twi.getColumnValueList('etayp')[istart]

        xtrack_twiss0 = _xtrack.TwissInit(betx = betax,
                                          alfx = alphax,
                                          dx = etax,
                                          dpx = etaxp,
                                          bety = betay,
                                          alfy = alphay,
                                          dy = etay,
                                          dpy = etayp)

        dict_twiss0 = {"betx":betax,
                       "alfx":alphax,
                       "dx":etax,
                       "dpx":etaxp,
                       "bety":betay,
                       "alfy":alphay,
                       "dy":etay,
                       "dpy":etayp}

        xtrack_line.set_particle_ref(pdg_id_0=11,
                                     p0c = m_e_eV*p0,
                                     s=s)

        xtrack_twiss = xtrack_line.twiss(method="4d",
                                     init=xtrack_twiss0)

    if elegant_ps is not None and elegant_twi is None:
        xtrack_particles = elegant2xsuite_particles(elegant_ps, xtrack_line)
        context = _xobjects.ContextCpu()
        xtrack_line.build_tracker(context)
        xtrack_line.track(xtrack_particles)

    return {"env":env,
            "xtrack_twiss0":xtrack_twiss0,
            "dict_twiss0":dict_twiss0,
            "xtrack_twiss":xtrack_twiss,
            "xtrack_particles":xtrack_particles}

def elegant2xsuite_particles(elegant_ps, xtrack_line) :
    if isinstance(elegant_ps, str):
        elegant_ps = _sdds.load(elegant_ps)

    x = _np.array(elegant_ps.getColumnValueList('x'))
    y = _np.array(elegant_ps.getColumnValueList('y'))
    xp = _np.array(elegant_ps.getColumnValueList('xp'))
    yp = _np.array(elegant_ps.getColumnValueList('yp'))
    t = _np.array(elegant_ps.getColumnValueList('t'))
    p = _np.array(elegant_ps.getColumnValueList('p'))


    xtrack_line.set_particle_ref(pdg_id_0=11,
                                 p0c=m_e_eV*p.mean())

    px = xp # TODO larger angle?
    py = yp # TODO larger angle?

    zeta = (t-t.mean())*_constants.c*xtrack_line.particle_ref.beta0[0]
    delta = (p-p.mean())/p.mean()

    particles = xtrack_line.build_particles(x=x,
                                            px=px,
                                            y=y,
                                            py=py,
                                            zeta=zeta,
                                            delta=delta)

    return particles

def xsuite2fbpic(particles) :
    '''Convert xsuite particles to fbpic dict of arrays'''
    pass

def fbpic2xsuite(fbpic_dict) :
    '''Convert fbpic dict of particle arrays to xsuite'''
    pass

def xsuite_Remove_DriftSlices(line) :
    '''Replace DriftSlices for Drifts and create deferred variable'''

    for ename in line.element_names :
        e = line[ename]
        if isinstance(e,_xtrack.DriftSlice) :
            length = e.weight * line[e.parent_name].length
            line.vars[ename+".L"] = length
            line.element_dict[ename] = _xtrack.Drift(length=line.vars[ename+".L"]._value)
            line.element_refs[ename].length = line.vars[ename + ".L"]


def xsuite_Add_ParticlesMonitor(line, num_particles=10000) :
    for ename in line.element_names :
        e = line[ename]

        pm = _xtrack.ParticlesMonitor(num_particles=num_particles,
                                      start_at_turn=0,
                                      stop_at_turn=1)
        line.insert(ename+".PM", pm, at=ename, anchor="end")

def xsuite_CalculateBeamSize(line) :

    beam = {}
    beam['element_names'] = []
    beam['n'] = []
    beam['s'] = []
    beam['x_mean'] = []
    beam['x_sigma'] = []
    beam['px_mean'] = []
    beam['px_sigma'] = []
    beam['y_mean'] = []
    beam['y_sigma'] = []
    beam['py_mean'] = []
    beam['py_sigma'] = []

    for ename in line.element_names :
        e = line[ename]

        if isinstance(e,_xtrack.ParticlesMonitor) :
            beam['element_names'].append(ename)
            beam['s'].append(e.s[0,0])
            beam['n'].append(len(e.x[:,0]))

            beam['x_mean'].append(e.x[:,0].mean())
            beam['x_sigma'].append(e.x[:,0].std())
            beam['px_mean'].append(e.px[:,0].mean())
            beam['px_sigma'].append(e.px[:,0].std())

            beam['y_mean'].append(e.y[:,0].mean())
            beam['y_sigma'].append(e.y[:,0].std())
            beam['py_mean'].append(e.py[:,0].mean())
            beam['py_sigma'].append(e.py[:,0].std())

    # make numpy arrays
    for k in beam.keys() :
        beam[k] = _np.array(beam[k])

    return beam
