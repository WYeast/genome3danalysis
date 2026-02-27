import numpy as np
import h5py
from . import _experimental_shape

AVAILABLE_SHAPES = ['sphere', 'ellipsoid', 'experimental']
    
def run(struct_id: int, hss_opt: h5py.File, params: dict) -> np.ndarray:
    """ Compute the lamin distance for a given structure.
    
    The lamin distance is the distance from each bead to the closest point in the lamina.
    
    The lamina is defined by a center and a shape: currently, only spheres and ellipsoids are supported.

    Args:
        struct_id (int): The index of the structure in the HSS file.
        hss_opt (h5py.File): The optimized HSS file, with coordinates of different structures in separate datasets.
        params (dict): A dictionary containing the parameters for the analysis.

    Returns:
        (np.ndarray): distances from each bead to the lamina.
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
        # Compute √(X² + Y² + Z²)
        rad_dist = np.linalg.norm(coord - center, axis=1)
        # Compute LamDist
        lamin_dist = rad_dist * (1 / rad_dist_scaled - 1)
        return lamin_dist
    # for experimental
    elif shape == 'experimental':
        return _experimental_shape.compute_lad_from_experimental_shape(struct_id, coord, params)
