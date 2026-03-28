import sdds as _sdds
from ._elegant import elegant_lte_loader as _elegant_lte_loader

try :
    import pybdsim
except ImportError :
    print("pybdsim not found, cannot convert elegant to bdsim line")

def elegant2bdsim_gmad(elegant_file,
                       line_name="FEBE",
                       start_element="CLA-FEA-MAG-QUAD-13",
                       end_element="CLA-FED-SIM-DUMP-01-START",
                       elegant_twi=None,
                       elegant_ps=None,
                       outputfilename="FEBE.gmad",
                       overwrite=True):

    # load elegant lattice file
    if isinstance(elegant_file, str) :
        lte = _elegant_lte_loader(elegant_file)

    # load elegant twiss file if provided
    if elegant_twi is not None :
        if isinstance(elegant_twi, str) :
            elegant_twi = _sdds.load(elegant_twi)

    # create bdsim machine
    machine = pybdsim.Builder.Machine()

    # loop over elements
    for k in lte :
        ee = lte[k] # elegant element

        # skip over LINE
        if ee['TYPE'].upper() == "LINE":
            continue

        if ee['TYPE'] == 'CHARGE' :
            machine.Append(pybdsim.Builder.Marker(k))
        elif ee['TYPE'] == 'DRIFT':
            machine.Append(pybdsim.Builder.Drift(k, l=ee['L']))

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
            subline_component_names.append(e)

        if e == end_element :
            adding = False

    machine.Write(outputfilename, overwrite=overwrite)

def elegant2bdsim_memory() :
    pass

def elegant2bdsim_particles() :
    pass