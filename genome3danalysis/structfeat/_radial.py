import numpy as np
import h5py
from . import _experimental_shape

AVAILABLE_SHAPES = ['sphere', 'ellipsoid', 'experimental']
AVAILABLE_EXP_RADIAL_MODES = ['max_lad', 'half_longest_axis']


def run(struct_id: int, h5_opt: h5py.File, params: dict) -> np.ndarray:
    """ Compute radial position of each bead in a structure.

    For sphere/ellipsoid this returns scaled distance to center.
    For experimental shape this returns a normalized radial-like score in [0, 1]:
        radial = 1 - clip(LAD / R, 0, 1)
    where R is chosen by ``exp_radial_mode``:
        - max_lad (default): R = max(LAD) per structure
        - half_longest_axis: R = 0.5 * longest axis of experimental shape
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
    coord = h5_opt['coordinates'][str(struct_id)][:]

    # compute radial distance
    # for sphere and ellipsoid
    if shape == 'sphere' or shape == 'ellipsoid':
        coord_scaled = np.array(coord) / np.array(radius)
        center_scaled = np.array(center) / np.array(radius)
        return np.linalg.norm(coord_scaled - center_scaled, axis=1)
    # for experimental
    elif shape == 'experimental':
        lad = _experimental_shape.compute_lad_from_experimental_shape(struct_id, coord, params)

        mode = params.get('exp_radial_mode', 'max_lad')
        assert mode in AVAILABLE_EXP_RADIAL_MODES, (
            f'exp_radial_mode must be one of {AVAILABLE_EXP_RADIAL_MODES}'
        )

        if mode == 'half_longest_axis':
            denom = _experimental_shape.compute_experimental_shape_half_longest_axis(struct_id, params)
        else:
            denom = np.max(lad)

        if denom <= 0:
            return np.ones_like(lad)

        norm_lad = np.clip(lad / denom, 0.0, 1.0)
        return 1 - norm_lad
