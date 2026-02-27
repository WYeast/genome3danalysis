import numpy as np
from alabtools.utils import Index
import h5py
from .. import utils

DEFAULT_DIST_CUTOFF = 500  # nm

def run(struct_id: int, hss_opt: h5py.File, params: dict) -> np.ndarray:
    """ Calculate the trans A/B ratio for each bead in the structure.
    
    Trans A/B ratio is defined as the ratio of the number of inter-chromosomal A beads within a distance threshold
    divided by the number of inter-chromosomal B beads within the same distance threshold:
        transAB[i] = N_trans_A[i] / (N_trans_A[i] + N_trans_B[i])
        where N_trans_A[i] is the number of inter-chromosomal A beads within a distance threshold of bead i,
        and N_trans_B[i] is the number of inter-chromosomal B beads within the same distance threshold of bead i.

    Args:
        struct_id (int): The index of the structure in the HSS file.
        hss_opt (h5py.File): The optimized HSS file, with coordinates of different structures in separate datasets.
        params (dict): A dictionary containing the parameters for the analysis.

    Returns:
        (np.ndarray): The trans A/B ratio for each bead in the structure.
    """
    
    # get the index
    index = Index(hss_opt)
    
    # Read A/B compartment file name
    try:
        filename = params['filename']
    except KeyError:
        raise KeyError('AB-compartment filename must be specified')
    # Load the A/B compartment file with pandas
    try:
        _, _, _, ab = utils.read_bed(filename, val_type=str)
    except KeyError:
        raise KeyError('AB-compartment file not found')
    # Assert that AB track is correct
    assert isinstance(ab, np.ndarray), 'AB-compartment track must be a numpy array'
    assert len(ab) == np.sum(index.copy == 0),\
        'AB-compartment track must have the same length as the number of haploid beads in the structure'
    assert len(np.unique(ab)) == 2 or len(np.unique(ab)) == 3,\
        'AB-compartment track must contain only A and B and optionally NA (unspecified format)'
    assert 'A' in np.unique(ab) and 'B' in np.unique(ab), 'AB-compartment track must contain A and B'    
    # Adapt AB track to multi-ploid index
    ab = utils.adapt_haploid_to_index(ab, index)
    
    # Get the surface-to-surface distance threshold
    try:
        dist_sts_thresh = params['dist_cutoff']
    except KeyError:
        dist_sts_thresh = DEFAULT_DIST_CUTOFF
    
    # get coordinates of struct_id
    coord = hss_opt['coordinates'][str(struct_id)][:]
    
    # get the radii of the beads
    radii = hss_opt['radii'][:]
    
    # Initialize the transAB_ratio
    transAB_ratio = np.zeros(len(index)).astype(float)
    
    # Loop over all beads
    for i in range(len(index)):
        
        # FIND PROXIMAL BEADS
        # First, we get the center-to-center distances between the bead i and all other beads
        dists = np.linalg.norm(coord - coord[i], axis=1)
        # Then, we get the center-to-center distance trhesholds between bead i and all other beads,
        # which is the sum of the radii and the surface-to-surface distance threshold:
        #       dcap_ij = ri + rj + d_sts_thresh, fixed i
        dcap = radii[i] + radii + dist_sts_thresh
        # Finally, we get the indices of the beads that are within the distance threshold
        prox_beads = np.where(dists < dcap)[0]
        # Remove the bead i from the proximal beads (no self-interactions)
        prox_beads = prox_beads[prox_beads != i]
        
        # FILTER INTER-CHROMOSOMAL FROM PROXIMAL BEADS
        # Get the chromosome and copy of the proximal beads
        chrom_prox_beads = index.chrom[prox_beads]
        copy_prox_beads = index.copy[prox_beads]
        # Get a mask that filters only the proximal beads that are inter-chromosomal (different chromosomes or different copies)
        inter_mask = np.logical_or(chrom_prox_beads != index.chrom[i], copy_prox_beads != index.copy[i])
        # Get the proximal beads that are inter-chromosomal
        prox_inter_beads = prox_beads[inter_mask]
        
        # GET TRANSAB RATIO
        # Get the A/B identity of the proximal inter-chromosomal beads
        ab_prox_inter_beads = ab[prox_inter_beads]
        # Get the TransAB ratio: n_trans(A) / (n_trans(A) + n_trans(B))
        n_trans_A = np.sum(ab_prox_inter_beads == 'A')
        n_trans_B = np.sum(ab_prox_inter_beads == 'B')
        transAB_ratio[i] = n_trans_A / (n_trans_A + n_trans_B)
        
        del dcap, dists, prox_beads, chrom_prox_beads, copy_prox_beads, inter_mask, prox_inter_beads, ab_prox_inter_beads
    
    return transAB_ratio
