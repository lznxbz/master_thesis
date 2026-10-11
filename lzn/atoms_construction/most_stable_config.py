def find_the_most_stable_config(in_dir, matching = None) -> list:
    from lzn.general import to_path
    from ase.io import read
    from ase.io.ulm import InvalidULMFileError

    in_dir = to_path(in_dir)

    # Find input files
    if isinstance(matching, str):
        in_list = sorted(in_dir.glob(matching))
    elif isinstance(matching, list):
        in_list = [in_dir / m for m in matching]
    elif matching is None:
        in_list = sorted(in_dir.glob("*.traj"))
    else:
        raise TypeError(f"Expected str or list, got {type(matching)}")

    out_list = []
    previous_site, previous_energy = None, None
    most_stable = None
    for in_file in in_list:
        site_name = in_file.stem.split("_")[-2]
        if site_name[:4] != "site":
            raise ValueError(f"Unrecognized filename: {in_file.name}. Expected '..._sitexx_cxx.traj'")
        try:
            atoms = read(in_file, index=-1, format="traj")
        except InvalidULMFileError:
            continue

        atoms_energy = atoms.get_potential_energy()
        if previous_site != site_name:
            out_list.append(most_stable)
            previous_site = site_name
            previous_energy = atoms_energy
            most_stable = in_file
        else:
            if atoms_energy < previous_energy:
                previous_energy = atoms_energy
                most_stable = in_file
            else:
                continue

    out_list.pop(0)
    return out_list
