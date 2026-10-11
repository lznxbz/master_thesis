from ase.io import read, write
from lzn.general import to_path, batch_do
from pathlib import Path



def asetraj_to_xyz(in_path : str | Path, out_dir : Path | str, center=False, comment=" ", *args):
    in_path = to_path(in_path)
    out_dir = to_path(out_dir)
    if in_path.suffix != ".traj":
        raise ValueError(f"File {in_path} is not a ASE trajactory file")
    
    file_atoms = read(in_path, format="traj") # Read the last frame by default
    if center:
        file_atoms.translate(-file_atoms.get_center_of_mass()) # center
    out_name = in_path.stem + ".xyz"
    out_path = out_dir / out_name
    if comment == "Auto":
        comment = f"Converted from {in_path}"
    write(out_path, file_atoms, format='xyz', comment=comment)

    return 0





if __name__ == "__main__":
    batch_do(
        asetraj_to_xyz,
        "~/TiO2_organized/771_r2scan_opt",
        "tio12CDa_*_opt-matpes_r2scan.traj",
        "~/TiO2_organized/781_r2scan_opt_xyz"
    )

