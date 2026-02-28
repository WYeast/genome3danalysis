import numpy as np
import h5py
from . import _experimental_shape

AVAILABLE_SHAPES = ['sphere', 'ellipsoid', 'experimental']
    
def run(struct_id: int, hss_opt: h5py.File, params: dict) -> np.ndarray:
    """ Compute the radial distance of each bead in a structure to a given center and shape.
    
    Currently, only sphere and ellipsoid shapes are supported.

    Args:
        struct_id (int): The index of the structure in the HSS file.
        hss_opt (h5py.File): The optimized HSS file, with coordinates of different structures in separate datasets.
        params (dict): A dictionary containing the parameters for the analysis.

    Returns:
        (np.ndarray): distances of each bead to the center of the shape.
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
    
    # compute radial distance
    # for sphere and ellipsoid
    if shape == 'sphere' or shape == 'ellipsoid':
        coord_scaled = np.array(coord) / np.array(radius)
        center_scaled = np.array(center) / np.array(radius)
        return np.linalg.norm(coord_scaled - center_scaled, axis=1)
    # for experimental
    elif shape == 'experimental':
        lad = _experimental_shape.compute_lad_from_experimental_shape(struct_id, coord, params)
        max_lad = np.max(lad)
        if max_lad <= 0:
            return np.ones_like(lad)
        norm_lad = lad / max_lad
        return 1 - norm_lad
