"""
To check whether pbe0_def2svp and omol25_low find the same local minimum.
By optimize (dft optimized structure) using mlip.
Then check RMSD minimum: if it close to 0
-> same local minimum.
-> the differences between two optimized structure caused by the method.
"""

from ase.io import read
from ase.optimize import BFGS

from sevenn.calculator import SevenNetCalculator

from lzn.general import *
from lzn.data_analysis.rmsd_min import find_rmsd_min





def compare_diff(
        in_path, # g16 log file
        out_dir, # save opt log and traj files
        ref_dir1,
        ref_dir2,
        calc,
    ):
    """

    """
    # Initiation
    in_path = to_path(in_path)
    out_dir = to_path(out_dir)
    ref_dir1 = to_path(ref_dir1)
    ref_dir2 = to_path(ref_dir2)

    # Step 0.5: parse filename, and get mlip opt name
    in_stem = in_path.stem
    structure_name = in_stem.split("_opt-")[0]
    othersite_number = int(structure_name[-2:]) + 1

    mlip_opt_suffix = "_opt-omol25_low.xyz"
    mlip_filename = structure_name + mlip_opt_suffix
    mlip_othersite_filename = "tio12CDa_co_site" + str(othersite_number) + mlip_opt_suffix

    # First step run: MLIP opt for g16 optimized structures
    atoms_compare = read(in_path, index=-1, format="gaussian-out")
    atoms_compare.calc = calc

    log_path = out_dir / (in_stem + ".log")
    traj_path = out_dir / (in_stem + ".traj")

    dyn = BFGS(atoms_compare, logfile=log_path, trajectory=traj_path, maxstep=0.2)
    dyn.run(fmax=0.02)


    # Second step load MLIP opt structures, and MLIP othersite opt structures
    mlip_filepath = ref_dir1 / mlip_filename
    mlip_othersite_filepath = ref_dir2 / mlip_othersite_filename

    atoms_ref = list(map(read, [mlip_filepath, mlip_othersite_filepath]))

    # Third step compare rmsd: load first, then cover by the next if rmsd is lower
    rmsd_info = list(map(find_rmsd_min, atoms_ref, [atoms_compare]*2))

    relavant_structure, rmsd_min = mlip_filename, rmsd_info[0][2]
    if rmsd_info[1][2] < rmsd_min: # change to while loop if more equivalent structures
        relavant_structure, rmsd_min = mlip_othersite_filename, rmsd_info[1][2]

    return relavant_structure, rmsd_min






if __name__ == "__main__":
    modal = "omol25_low"
    calc = SevenNetCalculator(model="7net-omni", modal=modal)
    out_list = batch_do(
        compare_diff,
        "~/TiO2_organized/07_catads_g16_opt",
        "tio12CDa_co_site*_opt-pbe0_def2svp.log",
        "~/TiO2_organized/09_catads_check_consistency",
        ref_dir1 = "~/TiO2_organized/08_catads_opt_xyz",
        ref_dir2 = "~/TiO2_organized/08_catads_opt_xyz",
        calc=calc,
    )
    print(out_list)

