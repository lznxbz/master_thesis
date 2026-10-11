from lzn.general import *
from ase.io import read


def extract_ase_energy(file_path: Path | str, unit="eV", *args) -> float:
    """
    Extract energy from an ASE traj file.
    Energy unit in hartree
        - file_path: Path to an ASE traj file
        - unit: eV (default), Hartree, kJ/mol
    """
    file_path = to_path(file_path)
    atoms = read(file_path, index= -1, format="traj")
    atoms_energy = atoms.get_potential_energy()

    match unit.lower():
        case "ev":
            return atoms_energy
        case "hartree":
            return atoms_energy * eV_to_hartree
        case "kj/mol":
            return atoms_energy * eV_to_kJdmol
        case _:
            raise ValueError(f"Unknown unit {unit}")



def extract_g16_energy(file_path: Path | str, unit="eV", *args) -> float:
    """
    Extract energy from a g16 log file.
    Energy unit in hartree
        - file_path: Path to a g16 log file
        - unit: eV (default), Hartree, kJ/mol
    """
    file_path = to_path(file_path)
    atoms = read(file_path, index= -1, format="gaussian-out")
    atoms_energy = atoms.get_potential_energy()

    match unit.lower():
        case "ev":
            return atoms_energy
        case "hartree":
            return atoms_energy * eV_to_hartree
        case "kj/mol":
            return atoms_energy * eV_to_kJdmol
        case _:
            raise ValueError(f"Unknown unit {unit}")




if __name__ == "__main__":
    ep = extract_g16_energy("~/TiO2_organized/07_catads_g16_opt/tio12CDa_co_site36_opt-pbe0_def2svp.log")
    print(ep)



