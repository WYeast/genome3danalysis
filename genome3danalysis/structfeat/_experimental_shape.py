import os
import struct
import numpy as np
from scipy import ndimage

try:
    import mrcfile
except ImportError:
    mrcfile = None


def _resolve_shape_file(struct_id: int, params: dict) -> str:
    if 'exp_shape_file' in params:
        shape_file = params['exp_shape_file'].format(struct_id=struct_id)
        if os.path.exists(shape_file):
            return shape_file

        # Backward-compatible fallbacks for nested per-structure directories,
        # e.g. /.../lamina_bin/{struct_id}/lamina/lamina.bin
        base, ext = os.path.splitext(shape_file)
        candidates = []
        if ext in ['.bin', '.mrc']:
            candidates.append(os.path.join(base, 'lamina', f'lamina{ext}'))
        else:
            candidates.extend([
                os.path.join(shape_file, 'lamina', 'lamina.bin'),
                os.path.join(shape_file, 'lamina', 'lamina.mrc')
            ])

        for cand in candidates:
            if os.path.exists(cand):
                return cand

        return shape_file

    if 'exp_shape_dir' in params:
        ext = params.get('exp_shape_ext', '.bin')
        default_file = os.path.join(params['exp_shape_dir'], f'{struct_id}{ext}')
        if os.path.exists(default_file):
            return default_file

        nested_file = os.path.join(params['exp_shape_dir'], str(struct_id), 'lamina', f'lamina{ext}')
        if os.path.exists(nested_file):
            return nested_file

        return default_file

    raise KeyError('Experimental shape requires exp_shape_file template or exp_shape_dir')


def _coords_to_voxel(coords: np.ndarray, origin: np.ndarray, dx: np.ndarray, shape_xyz: tuple) -> np.ndarray:
    vox = np.rint((coords - origin[None, :]) / dx[None, :]).astype(int)
    vox[:, 0] = np.clip(vox[:, 0], 0, shape_xyz[0] - 1)
    vox[:, 1] = np.clip(vox[:, 1], 0, shape_xyz[1] - 1)
    vox[:, 2] = np.clip(vox[:, 2], 0, shape_xyz[2] - 1)
    return vox


def _lad_from_bin(coords: np.ndarray, bin_name: str) -> np.ndarray:
    with open(bin_name, 'rb') as f:
        _, nx, ny, nz = struct.unpack('iiii', f.read(16))
        _ = struct.unpack('fff', f.read(12))  # center, currently unused
        origin = np.array(struct.unpack('fff', f.read(12)), dtype=float)
        dx = np.array(struct.unpack('fff', f.read(12)), dtype=float)
        raw = np.fromfile(f, dtype=np.int32)

    mappa = raw.reshape((nx, ny, nz, 4), order='C')
    vox = _coords_to_voxel(coords, origin, dx, (nx, ny, nz))

    bnd_idx = mappa[vox[:, 0], vox[:, 1], vox[:, 2], 0:3].astype(float)
    delta = (vox.astype(float) - bnd_idx) * dx[None, :]
    return np.linalg.norm(delta, axis=1)


def _lad_from_mrc(coords: np.ndarray, mrc_name: str, params: dict) -> np.ndarray:
    if mrcfile is None:
        raise ImportError('mrcfile is required to use experimental MRC nuclear shapes')

    with mrcfile.open(mrc_name, mode='r', permissive=True) as mrc:
        data = mrc.data
        transpose_mode = params.get('transpose_mode', 'none')
        if transpose_mode == 'T':
            data = data.T
        elif transpose_mode == 'swap_xz':
            data = np.swapaxes(data, 0, 2)

        origin = np.array([
            float(mrc.header.origin.x),
            float(mrc.header.origin.y),
            float(mrc.header.origin.z)
        ], dtype=float)

        if 'voxel_size' in params:
            dx = np.array(params['voxel_size'], dtype=float)
        else:
            dx = np.array([
                float(mrc.voxel_size.x),
                float(mrc.voxel_size.y),
                float(mrc.voxel_size.z)
            ], dtype=float)

    mask = (data > float(params.get('exp_shape_threshold', 0.5))).astype(np.uint8)
    edt_inside = ndimage.distance_transform_edt(mask, sampling=dx)
    edt_outside = ndimage.distance_transform_edt(1 - mask, sampling=dx)
    lad_map = np.where(mask > 0, edt_inside, edt_outside)

    vox = _coords_to_voxel(coords, origin, dx, mask.shape)
    return lad_map[vox[:, 0], vox[:, 1], vox[:, 2]]


def compute_lad_from_experimental_shape(struct_id: int, coords: np.ndarray, params: dict) -> np.ndarray:
    shape_file = _resolve_shape_file(struct_id, params)
    if shape_file.endswith('.bin'):
        return _lad_from_bin(coords, shape_file)
    if shape_file.endswith('.mrc'):
        return _lad_from_mrc(coords, shape_file, params)
    raise ValueError('Experimental shape file must be .bin or .mrc')
