from lzn.general import *

hardware_serial = f"%mem=16GB\n%nprocshared=8"
hardware_2parallel = f"%mem=8GB\n%nprocshared=4"

big_system_opt = "# opt=(gic,loose) freq pbe1pbe/def2svp empiricaldispersion=gd3bj"
small_system_opt = "# opt=(gic,tight) freq pbe1pbe/def2tzvp empiricaldispersion=gd3bj"

"""
Case
1. xyz -> g16 opt: head + copy coordinations, e.g. opt article structures, opt catads
2. opt chk -> g16 sp: head + read coord from chk, e.g. tzvp get more accurate energy for svp structures (clusters, catads)
3. sp chk -> g16 opt: head + read coordinations MO from chk, No such case for now
"""


### Read structure from g16 opt chk by g16, then do sp ###
gjf_g16opt_sp = lambda hardware, old_chk_path, chk_stem,comment: f"""{hardware}
%oldchk={old_chk_path}
%chk={chk_stem}.chk
# sp pbe1pbe/def2tzvp empiricaldispersion=gd3bj geom=allcheck guess=read

{comment}



"""

### Read mlip opt xyz file, and write to g16 sp gjf ###
gjf_mlipopt_sp = lambda hardware, chk_stem, comment, coordinations: f"""{hardware}
%chk={chk_stem}.chk
# sp pbe1pbe/def2tzvp empiricaldispersion=gd3bj

{comment}

0 1
{coordinations}


"""




def write_mlipopt_sp(
        from_file: Path | str,
        to_dir: Path | str,
        hardware: str,
        **kwargs
    ):
    """
    Read mlip opt xyz and write it to gjf. SAVE to harddrive !!!
        - from_file: XYZ file!!! mlip opt xyz
        - to_dir: g16 sp run folder
        - hardware: choose from serial or parallel
        - **kwargs: including all other parameters passed in, and not used.
    """
    from_file = to_path(from_file)
    if from_file.suffix != ".xyz":
        raise ValueError(f"use XYZ file!!! mlip opt xyz")

    from_stem = from_file.stem
    to_dir = to_path(to_dir)
    to_stem = from_stem + "_sp-pbe0_def2tzvp"
    to_name = to_stem + ".gjf"
    out_path = to_dir / to_name

    gjf_structure = ''.join(from_file.read_text(encoding="utf-8").splitlines(keepends=True)[2:])

    gjf_tot = gjf_mlipopt_sp(hardware, to_stem, f"from xyz file {from_file.name}", gjf_structure)

    out_path.write_text(gjf_tot, encoding="utf-8")
    print(f"Writen {out_path.name}")
    
    return out_path.name

def write_g16opt_sp(
        from_file: Path | str,
        to_dir: Path | str,
        hardware: str,
        **kwargs
    ):
    """
    Read g16 opt xyz and write it to gjf. SAVE to hard drive !!!
        - from_file: CHK file!!! g16 checkpoint files
        - to_dir: g16 run folder
        - hardware: choose from serial or parallel
        - **kwargs: includinng all other parameters passed in, and not used.
    """
    from_file = to_path(from_file)
    if from_file.suffix != ".chk":
        raise ValueError(f"use CHK file!!! g16 checkpoint file")
    from_stem = from_file.stem
    to_dir = to_path(to_dir)
    to_stem = from_stem + "_sp-pbe0_def2tzvp"
    to_name = to_stem + ".gjf"
    out_path = to_dir / to_name

    old_chk_path = str(from_file.absolute())
    gjf_tot = gjf_g16opt_sp(hardware, old_chk_path, to_stem, f"from chk file {from_file.name}" )

    out_path.write_text(gjf_tot, encoding="utf-8")
    print(f"Writen {out_path.name}")

    return out_path.name



def write_xyz_g16opt(
        from_file: Path | str,
        to_dir: Path | str,
        hardware: str,
        job: str,
        cm: str = "",
        **kwargs
    ):
    """
    Description required


    """
    from_file = to_path(from_file)
    if from_file.suffix != ".xyz":
        raise ValueError(f"File {from_file.name} does not end with .xyz")
    from_stem = from_file.stem
    to_dir = to_path(to_dir)

    # Naming of to_file
    name_dict = parse_lzn_stem(from_stem)
    name_dict["opt"] = "pbe0_def2svp"
    to_stem = f"{name_dict["structure_name"]}_opt-{name_dict["opt"]}"
    if cm:
        to_stem += "_cm-" + cm
    to_name = to_stem + ".gjf"
    out_path = to_dir / to_name

    gjf_structure = ''.join(from_file.read_text(encoding="utf-8").splitlines(keepends=True)[2:])
    gjf_tot = f"""{hardware}
%chk={to_stem}.chk
{job}

from xyz file {from_file.name}

0 1
{gjf_structure}


"""

    out_path.write_text(gjf_tot, encoding="utf-8")
    print(f"Writen {out_path.name}")

    return out_path.name




"""
Choose from:
write_mlipopt_sp, 
write_g16opt_sp, 
write_xyz_g16opt,
"""

if __name__ == "__main__":
    batch_do(
        write_mlipopt_sp,
        "~/TiO2_organized/04_adsorbates_mlip_opt",
        "*_r2scan.xyz",
        "/home/zining/TiO2_organized/041_adsorbates_sp",
        hardware = hardware_serial,
        job = big_system_opt,
    )



