import numpy as np
from scipy.spatial.distance import cdist
import h5py
from . import _experimental_shape

AVAILABLE_SHAPES = ['sphere', 'ellipsoid', 'experimental']
DEFAULT_EXTERIOR_THRESHOLD = 0.85
DEFAULT_TSA_EXPONENT = 0.004  # TODO: check this value
    
def run(struct_id: int, hss_opt: h5py.File, params: dict) -> np.ndarray:
    """ Compute the lamin TSA-seq signal for a given structure.
    
    The lamina is estimated from the most external beads of the structure.
    The TSA-seq signal is computed as the sum of the exponential of the distances to the lamina:
        TSA[i] = sum_j exp(-alpha * dist(i, j))

    Args:
        struct_id (int): The index of the structure in the HSS file.
        hss_opt (h5py.File): The optimized HSS file, with coordinates of different structures in separate datasets.
        params (dict): A dictionary containing the parameters for the analysis.

    Returns:
        (np.ndarray): TSA-distances to the lamina for each bead in the structure.
    """
    
    # Read the center or set it to [0, 0, 0]
    try:
        center = params['center']
    except KeyError:
        center = [0, 0, 0]
    assert len(center) == 3, 'Center must be a 3D vector'
    
    # Read the shape
    try:
        shape = params['shape']
    except KeyError:
        raise KeyError('Shape must be specified')
    assert shape in AVAILABLE_SHAPES, 'Shape not available'
    
    # Read the radius (only for sphere and ellipsoid)
    if shape == 'sphere' or shape == 'ellipsoid':
        try:
            radius = params['radius']
        except KeyError:
            raise KeyError('Radius must be specified')
        if shape == 'sphere':
            assert isinstance(radius, (int, float)), 'Radius must be a number since shape is a sphere'
        elif shape == 'ellipsoid':
            assert len(radius) == 3, 'Radius must be a 3D vector since shape is an ellipsoid'
            for r in radius:
                assert isinstance(r, (int, float)), 'Radius must be a 3D vector of numbers since shape is an ellipsoid'
    
    # Read the exterior threshold for lamina beads
    try:
        ext_thr = params['exterior_threshold']
    except KeyError:
        ext_thr = DEFAULT_EXTERIOR_THRESHOLD
    # Assert it's a float between 0 and 1
    assert isinstance(ext_thr, float), 'Exterior threshold must be a float'
    assert ext_thr >= 0 and ext_thr <= 1, 'Exterior threshold must be between 0 and 1'
    
    # Read the TSA-seq exponent
    try:
        tsa_alpha = params['tsa_exponent']
    except KeyError:
        tsa_alpha = DEFAULT_TSA_EXPONENT
    assert isinstance(tsa_alpha, float), 'TSA-seq exponent must be a float'
    
    # get coordinates of struct_id
    coord = hss_opt['coordinates'][str(struct_id)][:]
    
    # compute lamin distance
    # for sphere and ellipsoid
    if shape == 'sphere' or shape == 'ellipsoid':
        # We use the formula:
        #       LamDist = √(X² + Y² + Z²) * (1 / √((X/a)² + (Y/b)² + (Z/c)²) - 1)
        # Compute X/a, Y/b, Z/c
        coord_scaled = np.array(coord) / np.array(radius)
        center_scaled = np.array(center) / np.array(radius)
        # Compute √((X/a)² + (Y/b)² + (Z/c)²)
        rad_dist_scaled = np.linalg.norm(coord_scaled - center_scaled, axis=1)
        # Find the exterior beads and their coordinates
        lamina_beads = np.where(rad_dist_scaled > np.quantile(rad_dist_scaled, ext_thr))[0]
        lamina_coord = coord[lamina_beads, :]
        # Compute distance from every bead of the model to every bead of the lamina
        dist_to_lamina = cdist(coord, lamina_coord, 'euclidean')  # (nbead, nbead_lamina)  TODO: CHECK: ASLI TOOK THE MINIMUM OF THESE DISTANCES?
        # The TSA signal is then taken as a negative exponential of the distances previously computed
        lamina_tsa = np.sum(np.exp(- tsa_alpha * dist_to_lamina), axis=1)  # (nbead,)  TODO: CHECK SUM, MEAN OR MIN
        return lamina_tsa
    
    # for experimental
    elif shape == 'experimental':
        lamina_dist = _experimental_shape.compute_lad_from_experimental_shape(struct_id, coord, params)
        return np.exp(- tsa_alpha * lamina_dist)
