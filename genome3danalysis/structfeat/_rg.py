import numpy as np
import alabtools.geo
from alabtools.utils import Index
import h5py

DEFAULT_WINDOW_SIZE = 5  # TODO: should be given in Mb and converted to number of beads using the resolution
    
def run(struct_id: int, hss_opt: h5py.File, params: dict) -> np.ndarray:
    """ Compute the radius of gyration for each bead in the structure.

    Args:
        struct_id (int): The index of the structure in the HSS file.
        hss_opt (h5py.File): The optimized HSS file, with coordinates of different structures in separate datasets.
        params (dict): A dictionary containing the parameters for the analysis.

    Returns:
        (np.ndarray): radius of gyration for each bead in the structure.
    """
    
    # Read the widow size for the gyration radius calculation
    try:
        window = params['window']
    except KeyError:
        window = DEFAULT_WINDOW_SIZE
    # If the window is not an odd number, add 1 to it
    if window % 2 == 0:
        window += 1
    assert window >= 3, 'Window size should be at least 3'
    assert isinstance(window, int), 'Window size should be an integer'
    
    # get coordinates of struct_id
    coord = hss_opt['coordinates'][str(struct_id)][:]
    
    # get the radii of the beads
    radii = hss_opt['radii'][:]
    
    # get the index
    index = Index(hss_opt)
    
    # get the number of beads
    nbead = hss_opt.attrs['nbead']
    assert nbead == len(index)
    
    # initialize the gyration radius array
    gyr = []
    
    # loop over the beads and calculate the gyration radius
    for i in range(nbead):
        # get start and end of the window
        s = i - int((window - 1) / 2)
        e = i + int((window - 1) / 2)
        # check if the window is out of the structure
        if s < 0 or e > nbead - 1:
            gyr.append(np.nan)
            continue
        # check if the window is out of the chromosome
        if index.chrom[s] != index.chrom[e]:
            gyr.append(np.nan)
            continue
        # calculate the gyration radius
        gyr.append(alabtools.geo.RadiusOfGyration(coord[s:e, :], radii[s:e]))
    
    return np.array(gyr)
