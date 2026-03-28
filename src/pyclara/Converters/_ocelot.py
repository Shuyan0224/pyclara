import sdds as _sdds
from ._elegant import elegant_lte_loader as _elegant_lte_loader

try:
    import ocelot
    from ocelot import MagneticLattice
    from ocelot import Drift, Marker, SBend, Quadrupole, Sextupole
    from ocelot.adaptors.elegant_lattice_converter import ElegantLatticeConverter

except ImportError:
    print("ocelot not found, cannot convert yaml to ocelot line")

def elegant2ocelot(elegant_file,
                   start_element="CLA-FEA-SIM-DIP-04-END",
                   end_element="CLA-FED-SIM-DUMP-01-START",
                   elegant_twi = None,
                   elegant_ps = None) :

    # load elegant lattice file
    if isinstance(elegant_file, str) :
        lte = _elegant_lte_loader(elegant_file)

    # load elegant twiss file if provided
    if elegant_twi is not None :
        if isinstance(elegant_twi, str) :
            elegant_twi = _sdds.load(elegant_twi)

    # load elegant phase space file if provided
    if elegant_ps is not None :
        if isinstance(elegant_ps, str) :
            elegant_ps = _sdds.load(elegant_ps)

    oe_list = []
    # loop over elements
    for k in lte :
        ee = lte[k] # elegant element

        # skip over LINE
        if ee['TYPE'].upper() == "LINE":
            continue

        oe = None

        etype = ee['TYPE']
        ename = ee['NAME']

        if etype == 'CHARGE' :
            oe = Marker(eid=ename)
        elif etype == 'DRIFT':
            if ee['L'] == 0 :
                oe = Marker(eid = ename)
            else :
                oe  = Drift(eid = ename,
                            l = ee['L'])
        elif etype == 'CSRDRIFT':
            if ee['L'] == 0 :
                oe = Marker(eid = ename)
            else :
                oe  = Drift(eid = ename,
                            l = ee['L'])
        elif etype == 'LSCDRIFT':
            if ee['L'] == 0 :
                oe = Marker(eid = ename)
            else :
                oe  = Drift(eid = ename,
                            l = ee['L'])
        elif etype == "CSRCSBEND" or etype == "CSBEND":
            float(ee['L'])
            float(ee['ANGLE'])
            oe = SBend(eid = ename,
                       l = ee['L'],
                       angle = ee['ANGLE'])
        elif etype == 'KQUAD':
            oe = Quadrupole(eid = ename,
                            l = ee['L'],
                            k1 = ee['K1'])
        elif etype == "KSEXT":
            oe = Sextupole(eid = ename,
                           l = ee['L'],
                           k2 = ee['K2'])
        elif etype == "KICKER": # TODO implement correctly
            oe = Drift(eid = ename,
                       l = ee['L'])
        elif etype == "ECOL":
            oe = Drift(eid = ename,
                       l = ee['L'])
        elif etype == "MAXAMP" :
            oe = Marker(eid = ename)
        elif etype == 'WATCH':
            oe = Marker(eid = ename)
        elif etype == 'MONI' :
            oe = Drift(eid = ename,
                       l = ee['L'])
        elif etype == "RFCW" : # TODO implement correctly
            oe = Drift(eid=ename,
                       l = ee['L'])
        elif etype == "RFDF" : # TODO implement correctly
            oe = Drift(eid=ename,
                       l = ee['L'])
        else :
            print("element type ", etype, " not recognised, skipping")

        if oe is not None :
            oe_list.append(oe)

    # loop over elements selecting from sub range
    oe_list_used = []

    adding = False
    istart = -1

    for i, e in enumerate(oe_list) :
        if e.element.id == start_element :
            adding = True
            istart = i

        if adding :
            oe_list_used.append(e)

        if e == end_element :
            adding = False

    ml = MagneticLattice(oe_list_used)

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

        ocelot_twiss0 = ocelot.Twiss(beta_x = betax,
                                     alpha_x = alphax,
                                     Dx = etax,
                                     Dxp = etaxp,
                                     beta_y = betay,
                                     alpha_y = alphay,
                                     Dy = etay,
                                     Dyp = etayp)


    ocelot_twiss = ocelot.twiss(ml,ocelot_twiss0)

    return {"ocelot_lattice":ml,
            "ocelot_twiss0":ocelot_twiss0,
            "ocelot_twiss":ocelot_twiss}



