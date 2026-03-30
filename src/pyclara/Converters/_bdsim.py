import sdds as _sdds
from ._elegant import elegant_lte_loader as _elegant_lte_loader

try :
    import pybdsim
except ImportError :
    print("pybdsim not found, cannot convert elegant to bdsim line")

def elegant2bdsim_gmad(elegant_file,
                       line_name="FEBE",
                       start_element="CLA-FEA-SIM-DIP-04-END",
                       end_element="CLA-FED-SIM-DUMP-01-START",
                       elegant_twi=None,
                       elegant_ps=None,
                       outputfilename="FEBE.gmad",
                       overwrite=True,
                       flipMagnets = -1,
                       particleType="e-",
                       particleRestMass=0.511,
                       emitNX = 3,
                       emitNY = 3):

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

    start_element = start_element.replace("-","_")
    end_element = end_element.replace("-","_")



    # load elegant twiss file if provided
    if elegant_twi is not None :
        if isinstance(elegant_twi, str) :
            elegant_twi = _sdds.load(elegant_twi)

    be_list = []

    # loop over elements
    for k in lte :
        ee = lte[k] # elegant element
        ename = ee['NAME'].replace('-','_')
        etype = ee['TYPE']

        # skip over LINE
        if etype.upper() == "LINE":
            continue

        be = None
        if etype == 'CHARGE' :
            be = pybdsim.Builder.Marker(ename)
        elif etype == 'DRIFT':
            if ee['L'] != 0 :
                be = pybdsim.Builder.Drift(ename, l=ee['L'])
            else :
                be = pybdsim.Builder.Marker(ename)
        elif etype == 'CSRDRIFT':
            if ee['L'] == 0 :
                be = pybdsim.Builder.Marker(ename)
            else :
                be  = pybdsim.Builder.Drift(ename,
                                            l = ee['L'])
        elif etype == 'LSCDRIFT':
            if ee['L'] == 0 :
                be = pybdsim.Builder.Marker(ename)
            else :
                be  = pybdsim.Builder.Drift(ename,
                                            l = ee['L'])
        elif etype == "CSRCSBEND" or etype == "CSBEND":
            be = pybdsim.Builder.SBend(ename,
                                       ee['L'],
                                       ee['ANGLE'],
                                       magnetGeometryType="none",
                                       e1 = ee['E1'],
                                       e2 = ee['E2'],
                                       fint=ee['FINT'],
                                       fintx=ee['FINT'],
                                       hgap=ee['HGAP']) # TODO fix for short magnets
        elif etype == 'KQUAD':
            be = pybdsim.Builder.Quadrupole(ename,
                                            ee['L'],
                                            flipMagnets*ee['K1'])
        elif etype == "KSEXT":
            be = pybdsim.Builder.Sextupole(ename,
                                           ee['L'],
                                           flipMagnets*ee['K2'])
        elif etype == "KICKER": # TODO implement correctly
            be = pybdsim.Builder.TKicker(ename,
                                        hkick=0,
                                        vkick=0,
                                        l=ee['L'])
        elif etype == "ECOL":
            if ee['L'] !=  0 :
                be = pybdsim.Builder.Drift(ename,
                                           ee['L'])
            else :
                be = pybdsim.Builder.Marker(ename)
        elif etype == "MAXAMP" :
            be = pybdsim.Builder.Marker(name = ename)
        elif etype == 'WATCH':
            be = pybdsim.Builder.Marker(ename)
        elif etype == 'MONI' :
            be = pybdsim.Builder.Drift(ename,
                                       ee['L'])
        elif etype == "RFCW" : # TODO implement correctly
            be = pybdsim.Builder.Drift(ename,
                                       ee['L'])
        elif etype == "RFDF" : # TODO implement correctly
            be = pybdsim.Builder.Drift(ename,
                                       ee['L'])
        else :
            print("element type ", etype, " not recognised, skipping")

        if be is not None :
            be_list.append(be)

    # need to be uppercase
    line_component_names = [e.upper() for e in lte[line_name]['LINE']]

    if start_element is None or start_element == "":
        start_element = line_component_names[0]
    if end_element is None or end_element == "":
        end_element = line_component_names[-1]

    # select components within range
    subline_components = []

    adding = False
    istart = -1

    for i, e in enumerate(be_list) :
        if e.name == start_element:
            adding = True
            istart = i

        if adding :
            subline_components.append(e)

        if e.name == end_element:
            adding = False

    # create bdsim machine
    machine = pybdsim.Builder.Machine()
    for be in subline_components :
        machine.Append(be)

    # beam
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

        b  = pybdsim.Beam.Beam()
        b.SetParticleType(particleType)
        b.SetEnergy(p0*particleRestMass,"MeV")
        b.SetDistributionType("gausstwiss")
        b.SetBetaX(betax)
        b.SetBetaY(betay)
        b.SetAlphaX(alphax)
        b.SetAlphaY(alphay)
        b.SetDispX(etax)
        b.SetDispY(etay)
        b.SetDispXP(etaxp)
        b.SetDispYP(etayp)
        b.SetEmittanceNX(emitNX)
        b.SetEmittanceNY(emitNY)
        b.SetSigmaE(0.001)

        machine.AddBeam(b)

    machine.AddSampler('all')

    machine.Write(outputfilename, overwrite=overwrite)

def elegant2bdsim_memory() :
    pass

def elegant2bdsim_particles() :
    pass