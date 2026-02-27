import numpy as np
from scipy.spatial.distance import cdist
import pickle
import h5py
from . import _mcl

DEFAULT_TSA_EXPONENT = 0.004  # TODO: check this parameter!

def run(struct_id: int, hss_opt: h5py.File, params: dict, what_to_measure: str, body_type: str) -> np.ndarray:
    """ Calculate either the distance or the TSA-distance between each bead and the bodies of interest, e.g. speckles.
    
    The distance is calculated as the minimum distance between each bead and each body.
    
    The TSA-distance is calculated as the sum of the exponential of the negative distance between each bead and each body:
           TSA[i] = sum_j exp(-alpha * dist[i,j])

    Args:
        struct_id (int): The index of the structure in the HSS file.
        hss_opt (h5py.File): The optimized HSS file, with coordinates of different structures in separate datasets.
        params (dict): A dictionary containing the parameters for the analysis.
        what_to_measure (str): Either 'dist' or 'tsa', to calculate the distance or the TSA-distance, respectively.

    Returns:
        (np.ndarray): distance or TSA-distance for each bead.
    """
    
    # Obtain body centers from either a precomputed file or an on-the-fly MCL clustering.
    if params.get('use_mcl', False):
        bodies = _mcl.compute_mcl_cluster_coms(struct_id, hss_opt, params, body_type)
    else:
        # Read file name
        try:
            filename = params['filename']
        except KeyError:
            raise KeyError('Body filename must be specified (or set use_mcl=true with MCL params)')
        # Load the body file with pickle
        try:
            bodies = pickle.load(open(filename, 'rb'))
        except KeyError:
            raise KeyError('Body file not found')
        # Take the bodies of the structure
        try:
            bodies = bodies[struct_id]
        except KeyError:
            raise KeyError('Structure {} not found in body file'.format(struct_id))

    # Assert that bodies is a correct numpy array
    assert isinstance(bodies, np.ndarray), 'Bodies must be a numpy array'
    assert bodies.ndim == 2, 'Bodies must be a 2D array'
    assert bodies.shape[1] == 3, 'Bodies second dimension must be 3 (x,y,z)'
    
    # If there are no bodies, return an array of NaNs

    
    # get coordinates of struct_id
    coord = hss_opt['coordinates'][str(struct_id)][:]
    
    # If there are no bodies, return an array of NaNs
    if bodies.shape[0] == 0:
        return np.full(coord.shape[0], np.nan)
        
    # Compute the distance between each bead and each body
    dist_bodies = cdist(coord, bodies)  # shape: (n_beads, n_bodies)
    
    # Return the minimum distance for each bead
    if what_to_measure == 'dist':
        return np.min(dist_bodies, axis=1)
    elif what_to_measure == 'tsa':
        try:
            tsa_alpha = params['tsa_exponent']
        except KeyError:
            tsa_alpha = DEFAULT_TSA_EXPONENT
        return np.sum(np.exp(- tsa_alpha * dist_bodies), axis=1)
    