#!/usr/bin/env python3
import numpy as np
from ase import Atoms
from scipy.spatial.transform import Rotation as R

def rmsd(P, Q):
    """RMSD between two (N,3) arrays."""
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    return np.sqrt(np.mean(np.sum((P - Q)**2, axis=1)))

def best_fit_rotation_and_rmsd(P_ref, P_mov, weights=None):
    """
    Find the rotation that best aligns P_mov onto P_ref (minimum RMSD),
    assuming both are already centered to the same COM.

    Parameters
    ----------
    P_ref : (N,3) array
        Reference coordinates.
    P_mov : (N,3) array
        Moving coordinates (will be rotated).
    weights : (N,) array-like, optional
        Optional weights (e.g., atomic masses).

    Returns
    -------
    rot : scipy.spatial.transform.Rotation
        Optimal rotation (apply to P_mov).
    P_aligned : (N,3) array
        Rotated P_mov.
    rmsd_min : float
        RMSD after alignment.
    """
    P_ref = np.asarray(P_ref, float)
    P_mov = np.asarray(P_mov, float)

    # SciPy returns a Rotation that aligns the second argument onto the first.
    rot, rssd = R.align_vectors(P_ref, P_mov, weights=weights)  # rssd is a related fit metric.
    P_aligned = rot.apply(P_mov)
    rmsd_min = rmsd(P_ref, P_aligned)
    return rot, P_aligned, rmsd_min

def find_rmsd_min(
        a: Atoms,
        b: Atoms,
    ) -> tuple[np.ndarray, Atoms, float]:
    """
    Find the miminum RMSD between two atoms objects.
    To check if they are identical. (rmsd < 1E-6 are considered identical)
    Returns rot matrix, aligned b atoms, and rmsd min.
        - a: reference atoms
        - b: atoms that to be compared.
    """
    a.translate(-a.get_center_of_mass())
    b.translate(-b.get_center_of_mass())

    rot, b_aligned_pos, rmsd_min = best_fit_rotation_and_rmsd(a.positions, b.positions)

    b.positions[:] = b_aligned_pos
    print("RMSD (Å):", rmsd_min)
    return rot.as_matrix(), b, rmsd_min





def find_rmsd_min_partial(
        a,
        b,
        char_l,
):
    """
    Calc partial RMSD between two atoms objects, instead of for the whole particle.
    Returns rot matrix, aligned b atoms, and rmsd min.
    E.g. for TiO2 cluster + adsorbate, partial = TiO3 slab + adsorbate molecule.
        - a: reference atoms
        - b: atoms that to be compared.
        - char: a list of indices of adsorbate, NN Ti and O.
    """
    a = a[char_l]
    b = b[char_l]

    return find_rmsd_min(a, b)

