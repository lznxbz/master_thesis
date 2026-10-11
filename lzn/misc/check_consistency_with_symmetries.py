from lzn.general import *
from lzn.data_analysis.find_symmetry_equivalent_sites import *
from lzn.data_analysis.rmsd_min import find_rmsd_min
from lzn.data_analysis.extract_energies import *

import numpy as np
import pandas as pd

from ase.io import read



def center_inversion_read(atoms_path, swap_list):
    """
    This function read Atoms from atoms_path and center inversion,
    followed by remapping to center inverted atom indices (since find_rmsd_min is indices based)
    Return inverted Atoms object and its energy in unit of eV.
        - atoms_path: path to atoms file
        - swap_list: list of indices to swap
    """
    atoms = read(atoms_path)
    energy = atoms.get_potential_energy()

    #  Step 1: center inversion
    atoms_inv_pos = atoms.positions * -1.
    atoms.positions[:] = atoms_inv_pos

    # Step before 2: get indices of adsorbates
    order_adsorbates = list(range(36, len(atoms)))

    # Step 2: swap indices
    order = [None] * 36 + order_adsorbates # generate an index list. Empty in (TiO2)12, plus adsorbates
    for (i, j) in swap_list:
        order[i], order[j] = j, i

    atoms = atoms[order]
    return atoms, energy


def compare_even_and_inv_odd(even_path, odd_path, swap_list):
    """
    Calc energy diff and rmsd_min between even and inverted-odd indices Vsmax sites catads.
    return energy in unit of eV
        - even_path: path to even indices Atoms
        - odd_path: path to odd indices Atoms, passed to center_inversion_read
        - swap_list: list of indices to swap, passed to center_inversion_read
    """
    print(f"{even_path.name}, {odd_path.name}")
    atoms_even = read(even_path)
    energy_even = atoms_even.get_potential_energy()
    atoms_odd, energy_odd = center_inversion_read(odd_path, swap_list)
    rot, b, rmsd_min = find_rmsd_min(atoms_even, atoms_odd)
    ddE = energy_odd - energy_even
    return rmsd_min, ddE



def main_compare_even_and_inv_odd(ads, opt, where):
    """
    A batch_do version of compare_even_and_inv_odd
    Output energy unit = kJ/mol
    Read the code before run. Although it is save (does not modify anything on disk)
    """
    # get indices swap list
    atoms_symmtry_ref = read(
        to_path("~/TiO2_organized/05_mspa_clusters/tio12CDa_opt-pbe0_def2tzvpS_hsM.xyz"),
        index=1
    )
    dict_equiv = symmetry_equiv_sites(
        atoms_symmtry_ref,
        symm_fxn_i,
    )
    swap_list = dict_equiv["Ti"] + dict_equiv["O"]

    # load files
    file_dir = to_path(f"~/TiO2_organized/{where}")
    even_file_list = [f"tio12CDa_{ads}_site{i}_opt-{opt}.traj" for i in range(36,58,2)]
    odd_file_list = [f"tio12CDa_{ads}_site{i}_opt-{opt}.traj" for i in range(37,58,2)]
    name_list = [f"{ads}-site{i},{i+1}" for i in range(36,58,2)]

    # aquire data
    rmsd_even_odd = []
    ddE_list = []
    for i in range(len(even_file_list)):
        even_path = file_dir / even_file_list[i]
        odd_path = file_dir / odd_file_list[i]

        rmsd_min, ddE = compare_even_and_inv_odd(even_path, odd_path, swap_list)
        rmsd_even_odd.append(rmsd_min)
        ddE_list.append(ddE * eV_to_kJdmol)

    # Write df_omol_even_odd dataframe
    df_omol_even_odd = pd.DataFrame(zip(name_list, ddE_list, rmsd_even_odd), columns=["site #","ddE kJ/mol","rmsd min"])
    df_omol_even_odd.loc[:,"abs ddE"] = abs(df_omol_even_odd.loc[:,"ddE kJ/mol"])

    return df_omol_even_odd


def main_compare_same_site_diff_init(path1, path2, tol=1.,):
    """

    """
    # get indices swap list
    atoms_symmtry_ref = read(
        to_path("~/TiO2_organized/05_mspa_clusters/tio12CDa_opt-pbe0_def2tzvpS_hsM.xyz"),
        index=1
    )
    dict_equiv = symmetry_equiv_sites(
        atoms_symmtry_ref,
        symm_fxn_i,
    )
    swap_list = dict_equiv["Ti"] + dict_equiv["O"]

    # main
    path1, path2 = to_path(path1), to_path(path2)
    atoms1, atoms2 = read(path1), read(path2)

    ddE = atoms2.get_potential_energy() - atoms1.get_potential_energy()
    _rot, _b, rmsd_min = find_rmsd_min(atoms1, atoms2)
    inv = 0

    # rmsd > 1. => not the same site
    if rmsd_min > tol:
        atoms2_i, _e2 = center_inversion_read(path2, swap_list)
        _rot, _b, rmsd_min_inv = find_rmsd_min(atoms1, atoms2_i)
        # rmsd_min_inv < rmsd_min => inverted atoms2 is better matching site
        if rmsd_min_inv < rmsd_min:
            rmsd_min = rmsd_min_inv
            inv = 1

    return rmsd_min, ddE * eV_to_kJdmol, inv



if __name__ == "xxx":
    ls_tot = []
    for ads in ["co", "nh3", "h2o", "h2s"]:
        print(ads)
        ls_tot.append(main_compare_even_and_inv_odd(ads, "omol25_low", "07_catads_mlip_opt_from_tzvp"))

    df_tot = pd.concat(ls_tot, axis=0)
    print(df_tot)

if __name__ == "__main__":
    filename = "tio12CDa_nh3_site40_opt-pbe0_def2svp.log"
    dir_1 = to_path("~/TiO2_organized/07_catads_g16_opt")
    dir_2 = to_path("~/TiO2_organized/07_catads_g16_opt_from_tzvp")

    path_1 = dir_1 / "tio12CDa_nh3_site42_opt-pbe0_def2svp.log"
    path_2 = dir_2 / filename

    print(main_compare_same_site_diff_init(path_1, path_2))

