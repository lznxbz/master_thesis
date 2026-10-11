from lzn.atoms_construction.my_linear_algebra import *

# Read Tore MSPA hsEmsr .out file
from lzn.data_analysis.find_symmetry_equivalent_sites import *



def parse_vs_tolist(file_path, vs="max"):
    """
    parse Vsmax or Vsmin from .out, to a list of tuples
    (index, coords, Vs, at_element, at_index)
    Default: Vsmax
    """
    from lzn.general import to_path

    file_path = to_path(file_path)
    results = []
    start_flag = False

    with open(file_path, 'r') as f:
        for line in f:
            # Step 1: find start section
            if f"Information from Vs{vs} search" in line:
                start_flag = True
                continue

            if not start_flag:
                continue

            # Stop if section ends (optional safeguard)
            if not line.strip():
                break

            parts = line.split()

            # Skip non-data lines
            # Valid lines should start with an integer index
            if not parts or not parts[0].isdigit():
                continue

            try:
                # Example format:
                # 1 $ Ti 25  124.67   0.000020   16.36  -42.79
                #   2.708   -4.313    1.109   Ti 25  1.687   O 23  2.448   88.7

                idx = int(parts[0]) - 1

                Vs = float(parts[4])

                # coordinates
                x = float(parts[8])
                y = float(parts[9])
                z = float(parts[10])

                # nearest atom info
                at_element = parts[2]
                at_index = int(parts[3]) - 1

                results.append(
                    (idx, np.array([x, y, z]), Vs, at_element, at_index)
                )

            except (ValueError, IndexError):
                # skip malformed lines
                continue

    return results


def parse_standard_orientation_to_atoms(file_path):
    """
    parse standard orientation .out to ASE Atoms
    """
    from ase import Atoms

    symbols = []
    positions = []

    in_block = False

    with open(file_path, "r") as f:
        for line in f:
            s = line.strip()

            # start of coordinate block
            if "Standard Orientation in Angstrom" in line:
                in_block = True
                continue

            if not in_block:
                continue

            # end of coordinate block
            if set(s) == {"*"} and s:
                break

            parts = line.split()
            if len(parts) < 5:
                continue

            # expected: Element  Index  x  y  z
            # example: Ti 25  2.2014  -2.7899  0.5911
            symbol = parts[0]

            # skip malformed lines if first field is not an element-like token
            if not symbol[0].isalpha():
                continue

            try:
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
            except ValueError:
                continue

            symbols.append(symbol)
            positions.append((x, y, z))

    return Atoms(symbols=symbols, positions=positions)


def parse_vs_todict(file_path, vs="min"):
    """
    parse Vsmin or Vsmax from .out to a dict
    {at_atom_index: [(at_atoms_symbol, index, coords, Vs), ...], ..., "overwrite":[...]}
    Default: Vsmin
    """
    results = {}
    in_block = False
    overwrite = []

    with open(file_path, "r") as f:
        for line in f:
            # start block
            if f"Information from Vs{vs} search" in line:
                in_block = True
                continue

            if not in_block:
                continue

            if not line.strip():
                break

            parts = line.split()

            # skip non-data lines
            if not parts or not parts[0].isdigit():
                continue

            try:
                # --- basic indices ---
                nr = int(parts[0]) - 1           # running index (nr-1)
                atom_symbol = parts[2]           # e.g. "O"
                atom_index = int(parts[3]) - 1   # e.g. 36 → 35

                Vs = float(parts[4])             # potential

                # coordinates
                coord = np.array([
                    float(parts[8]),
                    float(parts[9]),
                    float(parts[10]),
                ])

            except (ValueError, IndexError):
                continue

            if atom_index not in results:
                results[atom_index] = (atom_symbol, nr, coord, Vs)
            else:
                results[atom_index] = [results[atom_index], (atom_symbol, nr, coord, Vs)]
                overwrite.append(atom_index)

    results["overwrite"] = overwrite

    return results


def parse_sym_equiv_sites(file_path, vs="max"):
    from lzn.general import to_path

    file_path = to_path(file_path)
    results = []
    start_flag = False

    with open(file_path, 'r') as f:
        sym_symbol0 = ""
        sym_equiv_sites = []
        for line in f:
            # Step 1: find start section
            if f"Information from Vs{vs} search" in line:
                start_flag = True
                continue
            if not start_flag:
                continue

            # Stop if section ends (optional safeguard)
            if not line.strip():
                results.append(tuple(sym_equiv_sites))
                break

            parts = line.split()
            # Skip non-data lines
            # Valid lines should start with an integer index
            if not parts or not parts[0].isdigit():
                continue
            try:
                # Example format:
                # 1 $ Ti 25  124.67   0.000020   16.36  -42.79
                #   2.708   -4.313    1.109   Ti 25  1.687   O 23  2.448   88.7
                idx = int(parts[0]) - 1
                sym_symbol1 = parts[1]
                if sym_symbol1 != sym_symbol0: # new site
                    if sym_symbol0 != "":
                        results.append(tuple(sym_equiv_sites))
                    sym_equiv_sites = [idx]
                else: # symmetry equivalent sites
                    sym_equiv_sites.append(idx)

                sym_symbol0 = sym_symbol1
            except (ValueError, IndexError):
                # skip malformed lines
                continue

    return {"X": results}



# Exclude symmetry equivalent sites (center inversion sites)
def exclude_unwanted_sites(file_path, find_sym="mspa.out"):
    """
    exclude center inversion sites and unstable sites (on oxygen)
    return vsmax_uni list of tuples, and (TiO2)n Atoms object.
    Kid of three predefined parse functions
    """
    from ase import Atoms
    # Step 1: Exclude center inversion sites
    atoms = parse_standard_orientation_to_atoms(file_path)

    vsmax = parse_vs_tolist(file_path)
    positions_vsmax = [i[1] for i in vsmax]
    x_vsmax = Atoms(f"X{len(positions_vsmax)}", positions_vsmax)

    atoms_ = atoms + x_vsmax

    if find_sym == "" or find_sym == "symm_fxn_i":
        dict_equiv = symmetry_equiv_sites(atoms_, symm_fxn_i)
        x_unique = [i[0] - len(atoms) for i in dict_equiv["X"]]
    elif find_sym == "mspa.out":
        dict_equiv = parse_sym_equiv_sites(file_path)
        x_unique = [i[0] for i in dict_equiv["X"]]
    else:
        raise ValueError(f"Unknown option find_sym: {find_sym}")

    vsmax_uni = [vsmax[i] for i in x_unique]
    # Step 2: Exclude unstable sites (on oxygen)
    vsmax_uni = [i for i in vsmax_uni if i[3] == "Ti"]

    return vsmax_uni, atoms


# Nearest O indices of each sites
def add_nearest_oxygen(filename):
    """
    Find three nearest oxygen atoms of input position.
    Kid of function exclude_unwanted_sites
    Return list of tuples: (x index, x coord, x val, NN type=Ti, NN index, nearest 3O index)
    and atoms
    """
    vsmax_uni_o = []

    vsmax_uni, atoms = exclude_unwanted_sites(filename)

    all_o_index = [i for i in range(len(atoms)) if atoms.symbols[i] == "O"]
    all_o_position = atoms.positions[all_o_index]

    for ii in vsmax_uni:
        pos = ii[1]

        all_o_dist = np.linalg.norm(all_o_position - pos, axis=1)
        fingerprint = np.argsort(all_o_dist)[:3]
        nearest_o_index = [all_o_index[i] for i in fingerprint]

        vsmax_uni_o.append(
            tuple(list(ii) + [nearest_o_index])
        )

    return vsmax_uni_o, atoms



# Branch: find the angle between O3 plane and Ti-Vsmax
def branch_find_angles_o3_tix(filename):
    """
    Branch kid of function add_nearest_oxygen
    Return list of angles(o3 plane, tix vec)
    """
    vsmax_uni_o, atoms = add_nearest_oxygen(filename)

    theta_list = []
    for ii in vsmax_uni_o:
        # Ti-Vsmax vector
        ti_coord = atoms.positions[ii[-2]]
        x_coord = ii[1]
        vec_tix = x_coord - ti_coord
        vec_tix = vec_tix / np.linalg.norm(vec_tix)

        # Normal vector of O3 plane
        o_positions = [atoms.positions[i] for i in ii[-1]]
        vo1 = o_positions[1] - o_positions[0]
        vo2 = o_positions[2] - o_positions[0]
        normal_vec = np.cross(vo1, vo2)
        normal_vec = normal_vec / np.linalg.norm(normal_vec)

        # Angle between two vectors
        theta = np.arccos(abs(np.dot(normal_vec, vec_tix)))
        theta_list.append(theta * 180 / np.pi)

    return theta_list



# Generate general binding configs
def calc_o3_plane_directions(filename) -> dict:
    """
    Generate six directions on o3 plane for featured sites, 60deg between each two.
    In the next step: align vector = main axis + o3_plane_direction
    Main axis could be normal vector of o3 plane OR TiX vector
    Kid of add_nearest_oxygen
    Return a dict:
    {
        "out":(x index, x coord, x val, NN type-Ti, NN index),
        "nno": [nearest 3 O indices],
        "dir_ino3":[six directions],
        "cat_vsmin":[Vsmins of dir_ino3 point toward],
    }
    and atoms
    """
    vsmax_, atoms = add_nearest_oxygen(filename)
    vsmin_ = parse_vs_todict(filename)

    vsmax_plus = {"out":[], "nno":[], "dir_ino3":[], "cat_vsmin":[]}

    for ii in vsmax_:
        o_positions = np.array([atoms.positions[i] for i in ii[-1]])

        # Get cluster Vsmin infos
        o_vsmin = []
        for i in ii[-1]:
            try:
                o_vsmin.append(vsmin_[i][-1])
            except(KeyError):
                o_vsmin.append(0.)

        # Get in O3 plane vectors (d00)
        o3_center = np.mean(o_positions, axis=0)

        o3_vecs = o_positions - o3_center
        o3_vecs_ = [normalize(i) for i in o3_vecs]
        o3_vecs = np.array(o3_vecs_)

        o3_6vecs = []
        o3_6vsmin = []
        for index_, i in enumerate(o3_vecs):
            o3_6vecs.append(i)
            o3_6vsmin.append(tuple([o_vsmin[index_]] * 2))
            if index_ == 2:
                o3_6vecs.append(normalize(i + o3_vecs[0]))
                o3_6vsmin.append(tuple([o_vsmin[index_], o_vsmin[0]]))
                break
            o3_6vecs.append(normalize(i + o3_vecs[index_ + 1]))
            o3_6vsmin.append(tuple([o_vsmin[index_], o_vsmin[index_ + 1]]))

        l_ = list(ii)
        oi_ = tuple(l_.pop(-1))
        vsmax_plus["out"].append(l_)
        vsmax_plus["nno"].append(oi_)
        vsmax_plus["dir_ino3"].append(o3_6vecs)
        vsmax_plus["cat_vsmin"].append(o3_6vsmin)

    return vsmax_plus, atoms



if __name__ == "__main__":
    # import pickle
    # filename = "/home/zining/organized/05_mspa_clusters/tio12CDa_opt-pbe0_def2tzvpS_hs.out"
    # list_info, _ = add_nearest_oxygen(filename)
    # print(len(list_info))
    # with open(
    #         "/home/zining/TiO2_organized/05_mspa_clusters/tio12CDa_opt-pbe0_def2tzvp_parse.pkl",
    #         "wb"
    # ) as f:
    #     pickle.dump(list_info, f)

    print(exclude_unwanted_sites(
        "/home/zining/organized/05_mspa_clusters/tio12CDa_opt-pbe0_def2tzvpS_hs.out",
        "symm_fxn_i"))
    print(exclude_unwanted_sites("/home/zining/organized/05_mspa_clusters/tio12CDa_opt-pbe0_def2tzvpS_hs.out"))

