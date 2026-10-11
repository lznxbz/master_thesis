#!/usr/bin/env python

from ase.io import read, write
from ase.optimize import BFGS

from sevenn.calculator import SevenNetCalculator

from lzn.general import to_path, parse_lzn_stem, write_lzn_stem, batch_do
import time


def overwrite_shift_of_calculator(calc, modal):
    model = calc.model
    model_state_dict = model.state_dict()
    shift_key = 'rescale_atomic_energy.shift'
    shift = model_state_dict[shift_key]

    config = calc.sevennet_config
    modal_map = config['_modal_map']
    zeros = torch.zeros(
        shift[modal_map[modal]].shape,
        dtype=shift[modal_map[modal]].dtype
    )

    # overwrite shift value
    shift[modal_map[modal]] = zeros
    # save it to state dict
    model_state_dict[shift_key] = shift
    # load model from overwritten state dict
    calc.model.load_state_dict(model_state_dict)

    return calc

# Preparation: read input files to Atoms objects, and naming output files
def filename_preparations_each(in_xyz, out_run_dir, out_xyz_dir, output_stem = "", calc_method=""):
    in_xyz = to_path(in_xyz)
    out_run_dir = to_path(out_run_dir)
    if out_xyz_dir:
        out_xyz_dir = to_path(out_xyz_dir)
    else:
        out_xyz_dir = out_run_dir

    # Output files
    if output_stem == "":
        dict_lzn = parse_lzn_stem(in_xyz) # parse information from filename to a dict
        dict_lzn["opt"] = calc_method
        stem = write_lzn_stem(dict_lzn)
    else:
        stem = output_stem

    log_path = out_run_dir / f"{stem}.log"
    traj_path = out_run_dir / f"{stem}.traj"
    out_xyz = out_xyz_dir / f"{stem}.xyz"

    return log_path, traj_path, out_xyz





def run_mlip_opt(
        in_xyz,
        out_run_dir,
        calc: SevenNetCalculator,
        out_xyz_dir = "", # out_run_dir by default
        fmax: float = 0.01,
        maxstep: float = 0.2,
        output_method: str = "save",
        output_stem: str = "",
        write_traj: bool = True,
    ):
    """
    run mlip optimization for a single structure, return initial and final energy (in eV)
    it will save log and traj files in out_run_dir, and save optimized xyz file in out_xyz_dir with comment containing energy info
    with filename in LZNxyz standard format
     - in_xyz: input structure (xyz file)
     - out_run_dir: directory to save log and traj files
     - out_xyz_dir: directory to save optimized xyz file, out_run_dir by default
     - calc: SevenNetCalculator object used for calculation
     - fmax: termination criteria - maximum force
     - maxstep: termination criteria - maximum step
     - output_method: options ("save" or "return") to save xyz file or return Atoms object
    """

    time_start = time.time()

    in_xyz = to_path(in_xyz)

    calc_method = calc.modal
    log_path, traj_path, out_xyz = filename_preparations_each(
        in_xyz, out_run_dir, out_xyz_dir, output_stem, calc_method)

    with in_xyz.open() as f:
        comment = f.readlines()[1].strip()


    """run optimization for a single structure with calc"""
    atoms = read(in_xyz)
    atoms.center(vacuum=10.)
    atoms.calc = calc

    traj_path = traj_path if write_traj else None # Save traj file or not

    e_initial = atoms.get_potential_energy()
    dyn = BFGS(atoms, logfile=log_path, trajectory=traj_path, maxstep=maxstep)
    dyn.run(fmax=fmax)
    e_final = atoms.get_potential_energy()

    atoms.translate(-atoms.get_center_of_mass()) # centerize

    print(f"finish {in_xyz.stem}, runtime = {time.time() - time_start}s")

    if output_method == "save":
        write(out_xyz, atoms, format='xyz', comment=comment)
        print(f"Save xyz file to {out_xyz}")
        return (e_initial, e_final) # unit is a.u.
    elif output_method == "return":
        print(f"finish opt, and return atoms object of {out_xyz.stem}")
        return atoms
    else:
        raise ValueError(f"{output_method} is not supported")


shift_calc = lambda calc: overwrite_shift_of_calculator(calc, calc.modal) if calc.modal == "omol25_low" else calc

def mlip_opt_batch(modal=""):
    time_init = time.time()
    in_dir = "/home/zining/TiO2_organized/08_catads_opt_xyz_from_tzvp"
    matching = "tio12CDa_*_opt-omol25_low.xyz"
    out_dir = "/home/zining/TiO2_organized/07_catads_mlip_opt_from_tzvp"
    out_xyz_dir="/home/zining/TiO2_organized/08_catads_opt_xyz_from_tzvp"

    calc = shift_calc(SevenNetCalculator(model="7net-omni", modal="matpes_r2scan"))

    batch_do(
        run_mlip_opt,
        in_dir,
        matching,
        out_dir,
        calc = calc,
        out_xyz_dir=out_xyz_dir,
        fmax=0.01,
    )

    print("DONE \\(-v-)/")
    print(f"Total time = {time.time() - time_init}")
    return 0

def mlip_opt_single():
    modal = "matpes_r2scan"
    calc = shift_calc(SevenNetCalculator(model="7net-omni", modal=modal))

    ads_l = ["co", "h2o", "h2s", "nh3"]

    for f in ["h2s", "nh3"]:
        run_mlip_opt(
            in_xyz=f"~/TiO2_organized/04_adsorbates_mlip_opt/{f}_opt-omol25_low.xyz",
            out_run_dir="~/TiO2_organized/04_adsorbates_mlip_opt",
            out_xyz_dir="~/TiO2_organized/04_adsorbates_mlip_opt",
            calc=calc,
            fmax=0.01,
            output_stem="",
        )
    return 0



if __name__ == "__main__":
    main_batch()

