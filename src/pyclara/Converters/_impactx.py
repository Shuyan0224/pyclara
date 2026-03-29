import sdds as _sdds
from ._elegant import elegant_lte_loader as _elegant_lte_loader

try :
    from impactx import ImpactX as _ImpactX
    from impactx import elements as _elements
    from impactx import twiss as _twiss
    from impactx import distribution as _distribution
except ImportError:
    print("impactx not found, cannot convert elegant to impactx line")

def elegant2impactx(elegant_file,
                    start_element="CLA-FEA-SIM-DIP-04-END",
                    end_element="CLA-FED-SIM-DUMP-01-START",
                    elegant_twi=None,
                    elegant_ps=None):
    '''
    :param elegant_file:
    :type elegant_file: str or dict
    :param start_element: start element (inclusive)
    :type start element: str
    :param end_element: end element (inclusive)
    :type end_element: str
    :param elegant_twi: Twiss SDDS file from Elegant
    :type elegant_twi: str or sdds elegant twi file
    :param elegant_ps: Watch SDDS file from Elegant
    :type elegant_ps: str or sdds elegant watch output file
    :return: dict of impactx objects
    '''
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

    # impactx elements
    ie_list = []

    # loop over elements
    for k in lte :
        ee = lte[k] # elegant element
        ename = ee['NAME']
        etype = ee['TYPE']

        # skip over LINE
        if etype.upper() == "LINE":
            continue

        ie = None

        if etype == 'CHARGE' :
            ie = _elements.Marker(name=ename)
        elif etype == "DRIFT":
            ie = _elements.Drift(name=ename, ds=ee['L'])
        elif etype == 'CSRDRIFT':
            if ee['L'] == 0 :
                ie = _elements.Marker(name = ename)
            else :
                ie  = _elements.Drift(name = ename,
                                      ds = ee['L'])
        elif etype == 'LSCDRIFT':
            if ee['L'] == 0 :
                ie = _elements.Marker(eid = ename)
            else :
                ie  = _elements.Drift(name = ename,
                                      ds = ee['L'])
        elif etype == "CSRCSBEND" or etype == "CSBEND":
            ie = _elements.ExactSbend(name = ename,
                                      ds = ee['L'],
                                      phi = ee['ANGLE'])
        elif etype == 'KQUAD':
            ie = _elements.ExactQuad(name = ename,
                                     ds = ee['L'],
                                     k = ee['K1'])
        elif etype == "KSEXT":
            ie = _elements.ExactMultipole(name = ename,
                                     ds = ee['L'],
                                     k_normal = [0,0,ee['K2']],
                                     k_skew = [0,0,0])
        elif etype == "KICKER": # TODO implement correctly
            ie = _elements.Drift(name = ename,
                                 ds = ee['L'])
        elif etype == "ECOL":
            ie = _elements.Drift(name = ename,
                                 ds = ee['L'])
        elif etype == "MAXAMP" :
            ie = _elements.Marker(name = ename)
        elif etype == 'WATCH':
            ie = _elements.Marker(name = ename)
        elif etype == 'MONI' :
            ie = _elements.Drift(name = ename,
                                 ds = ee['L'])
        elif etype == "RFCW" : # TODO implement correctly
            ie = _elements.Drift(name=ename,
                                 ds = ee['L'])
        elif etype == "RFDF" : # TODO implement correctly
            ie = _elements.Drift(name=ename,
                                 ds = ee['L'])
        else :
            print("element type ", etype, " not recognised, skipping")

        if ie is not None :
            ie_list.append(ie)

    # loop over elements selecting from sub range
    ie_list_used = []

    adding = False
    istart = -1

    for i, e in enumerate(ie_list) :
        if e.name == start_element :
            adding = True
            istart = i

        if adding :
            ie_list_used.append(e)

        if e.name == end_element :
            adding = False

    # get charge
    bunch_charge_C = lte['START']['TOTAL']
    npart = 50000

    impactx_twiss0 = None
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

        impactx_twiss0 = _twiss(beta_x=betax,
                                beta_y=betay,
                                beta_t=0.5,
                                emitt_x=3e-6,
                                emitt_y=3e-6,
                                emitt_t=2e-06,
                                alpha_x=alphax,
                                alpha_y=alphay,
                                alpha_t=0.0,
                                dispersion_x=etax,
                                dispersion_y=etay,
                                dispersion_px=etaxp,
                                dispersion_py=etayp,)

    # require a twiss file
    if elegant_twi is None :
        raise ValueError("elegant twiss file not provided, cannot set initial distribution for impactx")

    # impactx simulation
    sim = _ImpactX()

    # set numerical parameters and IO control
    sim.space_charge = False
    # sim.diagnostics = False  # benchmarking
    sim.slice_step_diagnostics = True

    # domain decomposition & space charge mesh
    sim.init_grids()

    #   reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(-1.0)
    ref.set_mass_MeV(0.510998950)
    ref.set_kin_energy_MeV(250)

    distr = _distribution.Waterbag(**impactx_twiss0)
    sim.add_particles(bunch_charge_C, distr, npart)

    # assign a fodo segment
    sim.lattice.extend(ie_list_used)

    # run simulation
    sim.track_particles()

    return {"impactx_sim":sim}
