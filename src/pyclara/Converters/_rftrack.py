import sdds as _sdds
from ._elegant import elegant_lte_loader as _elegant_lte_loader

try:
    import RF_Track
except ImportError:
    print("RF_Track not found, cannot convert elegant to RF_Track line")

def elegant2rftrack(elegant_file,
                    line_name="FEBE",
                    start_element="CLA-FEA-MAG-QUAD-13",
                    end_element="CLA-FED-SIM-DUMP-01-START",
                    elegant_twi=None):

    # load elegant lattice file
    if isinstance(elegant_file, str) :
        lte = _elegant_lte_loader(elegant_file)

