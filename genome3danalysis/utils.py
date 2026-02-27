import os
import numpy as np
import pandas as pd
import h5py
from alabtools.utils import Index
from alabtools.analysis import HssFile


def read_bed(filename: str, val_type: type = float) -> tuple:
    """ Read a bed file and return the chroms, starts, ends, and vals.
    
    The type of the fourth column can be specified with val_type.
    
    If the bed file does not contain a fourth column, chrom/starts/ends are returned.

    Args:
        filename (str): name (and path) of the bed file to read.
        val_type (type, optional): type of the fourth column. Defaults to float.

    Returns:
        (np.ndarray): chroms, of type str.
        (np.ndarray): starts, of type int.
        (np.ndarray): ends, of type int.
        (np.ndarray, optional): vals, of type val_type.
    """
    # Read the bed file using pandas
    # The separation key '\s+' means that the file is separated by any number of spaces
    # The header is set to None, so that the first line is not considered as the header
    bed = pd.read_csv(filename, sep='\s+', header=None)
    # Unpack chroms, starts, ends with the correct types
    chroms = bed[0].values.astype(str)
    starts = bed[1].values.astype(int)
    ends = bed[2].values.astype(int)
    # Try to unpack vals with the correct type
    try:
        vals = bed[3].values.astype(val_type)
    # If the file does not contain a fourth column, return chroms, starts, ends
    except IndexError:
        return chroms, starts, ends
    # Otherwise, return chroms, starts, ends, vals
    return chroms, starts, ends, vals


def adapt_haploid_to_index(hap_track: np.ndarray, index: Index) -> np.ndarray:
    """Given a haploid track, adapts it to the multi-ploid index.

    Args:
        hap_track (np.ndarray(nbead_hap)): haplotype track
        index (alabtools.utils.Index): index object
    
    Return:
        np.ndarray(nbeads): adapted track
    """
    multi_track = np.zeros(len(index), dtype=hap_track.dtype)
    copy_index = index.copy_index
    for i in copy_index:
        # copy_index is a Dictionary, where:
        # - keys are haploid indices (0, 1, 2, ..., nbead_hap - 1)
        # - values are the multiploid indices for the corresponding haploid index
        # for examples {0: [0, 1000], 1: [1, 1001], ...}
        multi_track[copy_index[i]] = hap_track[i]
    return multi_track

def parse_gap_track(filename: str) -> np.ndarray:
    """Read a 4-column BED gap/domain track and return a boolean gap mask.

    Supported 4th-column formats:
      - booleans / boolean-like strings: True/False, 1/0, yes/no, y/n
      - domain labels: regions labeled as domain/dom are treated as non-gap,
        all other labels are treated as gaps (e.g. cen/tel/gap/contig).

    Args:
        filename (str): BED file path.

    Returns:
        np.ndarray(bool): Haploid gap mask where True means gap/non-domain.
    """

    _, _, _, vals = read_bed(filename, val_type=str)
    vals = np.asarray(vals).astype(str)

    true_tokens = {'true', '1', 'yes', 'y', 't'}
    false_tokens = {'false', '0', 'no', 'n', 'f'}
    domain_tokens = {'domain', 'dom'}

    gap_mask = np.zeros(len(vals), dtype=bool)
    for i, raw in enumerate(vals):
        tok = raw.strip().lower()
        if tok in true_tokens:
            gap_mask[i] = True
        elif tok in false_tokens:
            gap_mask[i] = False
        elif tok in domain_tokens:
            gap_mask[i] = False
        else:
            gap_mask[i] = True

    return gap_mask


def create_optimized_hss(out_name: str, hss: HssFile, overwrite: bool = False) -> None:
    """ Creates an optimized HSS file, where the data of each structure is stored as a separate dataset.
    
    The file contains the following datasets/groups:
        - index (group): a group containing the index of the HSS file
        - coordinates (group): a group containing the coordinates of the structures
            *) each structure ('0', '1', '2', ...) is stored as a dataset
        - radii (dataset): the radii of the structures
    and the following attributes:
        - nbead (int): number of beads
        - nstruct (int): number of structures

    Args:
        hss (HssFile): HSS file to optimize
        overwrite (bool, optional): whether to overwrite an existing output file.
        Defaults to False.
    """
    
    # Make sure that the out_name has a valid path
    if not os.path.exists(os.path.dirname(out_name)):
        raise ValueError(f"Path {os.path.dirname(out_name)} does not exist.")
    # Make sure that the file does not already exist unless overwrite=True
    if os.path.exists(out_name):
        if overwrite:
            os.remove(out_name)
        else:
            raise ValueError(
                f"File {out_name} already exists. "
                "Delete it or call create_optimized_hss(..., overwrite=True)."
            )
    
    # Create a the HDF5 file for the optimized HSS
    h5 = h5py.File(out_name, 'w')
    
    # Add the attributes
    h5.attrs['nbead'] = hss.nbead
    h5.attrs['nstruct'] = hss.nstruct
    
    # Save the index
    hss.index.save(h5)
    
    # Create a dataset for the radii
    h5.create_dataset('radii', data=hss.radii, dtype=hss.radii.dtype)
    
    # Load the coordinates to memory from the HSS file
    coord = hss.coordinates
    
    # Create a group for the coordinates
    h5.create_group('coordinates')
    
    # Loop over the structures and save the coordinates
    for structID in range(hss.nstruct):
        h5['coordinates'].create_dataset(str(structID), data=coord[:, structID, :], dtype=coord.dtype)
    
    # Close the file
    h5.close()
