import numpy as np

calc_distance = lambda a,b: np.sum((a - b)**2)


def symm_fxn_i(
        indices_atom: list[int],
        positions_atoms: np.ndarray,
        center: np.ndarray,
        tol: float = 1e-3,
        mute = True
    ) -> list[tuple]:
    """
    Get inversion equivalent sites: list of tuples (index1,index2,...)
    print info, if distance is greater than tol
        - indices_atom: list of indices of atoms
        - positions_atoms: positions of atoms
        - center: center of symmetry
        - tol: float tolerance, default 1e-3
    """
    if __name__ == "__main__":
        mute = False

    # Centerize to inversion center
    positions_atoms -= center

    # Find equiv sites
    equiv_list = []
    while len(indices_atom) > 1:
        operated_atom_index = indices_atom.pop(0)
        operated_atom_coords = - positions_atoms[operated_atom_index] # apply inversion operation

        # Find the atom with min distance to operated_atom_coords
        min_distance = None
        min_index = None
        for compare_atom_index in indices_atom:
            current_distance = calc_distance(positions_atoms[compare_atom_index], operated_atom_coords)
            if (min_distance == None):
                min_distance, min_index = current_distance, compare_atom_index
                continue
            elif (current_distance < min_distance):
                min_distance, min_index = current_distance, compare_atom_index

        equiv_list.append((operated_atom_index, min_index))
        indices_atom.remove(min_index)

        # Check tolerance
        if (min_distance > tol) & (not mute):
            print(f"Exceed tolerance, distance = {min_distance}, between atom {operated_atom_index} and {min_index}")

    return equiv_list




def symmetry_equiv_sites(
        atoms,
        symm_fxn: callable,
        tol: float = 1e-3,
        center_of_symmetry: np.ndarray | str = np.array([0., 0., 0.]),
    ) -> dict:
    """
    Recieves an Atoms object with given symmetry,
    returns a dict {atom1:[(index1,index2,...),...], atom2:[...], ...},
    where tuples (index1,index2,...) record the indices of equivalent atoms under given symmetry.
        - atoms: ase.Atoms object with given symmetry
        - symm_fxn: callable symmetry finding function
        - tol: float tolerance, default 1e-3
        - center_of_symmetry: np.ndarray or str, default np.array([0., 0., 0.])
    """
    # centerize
    if isinstance(center_of_symmetry, np.ndarray):
        center = center_of_symmetry
    elif center_of_symmetry in [None, "", "center of mass"]:
        center = atoms.get_center_of_mass()
    else:
        raise ValueError(f"Unrecognized center_of_symmetry: {center_of_symmetry}. Input: np.ndarray or str")

    set_atoms = set(atoms.get_chemical_symbols())
    dict_equiv = dict(zip(set_atoms, [[]]*len(set_atoms)))
    positions_atoms = atoms.get_positions()

    # separate atom types
    for atom in set_atoms:
        indices_atom = [i for i, _atom in enumerate(atoms) if _atom.symbol == atom]
        dict_equiv[atom] = symm_fxn(indices_atom, positions_atoms, center, tol)

    return dict_equiv



if __name__ == "__main__":
    from ase.io import read
    from lzn.general import to_path

    atoms = read(
        to_path("~/organized/05_mspa_clusters/tio12CDa_opt-pbe0_def2tzvpS_hsM.xyz"),
        index=1,
    )

    dict_equiv = symmetry_equiv_sites(
        atoms,
        symm_fxn_i,
        tol=1e-3,
        center_of_symmetry=np.array([0., 0., 0.]),
    )

    print(dict_equiv)



