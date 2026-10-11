from lzn.general import to_path, parse_lzn_stem
from lzn.atoms_construction.parse_hs95_out import parse_standard_orientation_to_atoms, parse_vs_tolist, calc_o3_plane_directions
import lzn.atoms_construction.known_ligands as known_ligands



"""
This version combined the benefits from both v1 and v2.
1. Align adsorbate with cluster using: Vsmin-Donor // Acceptor-Vsmax, which is the same as version 1.
2. Apply tilt on adsorbates, which is introduced in version 2.

This code used for build catalyst-adsorbates initiation structure for (TiO2)n clusters on the featured sites. 
featured sites: Vsmax at the center of TiO3 plate, atop configuration on Ti.

Try 6 possible binding configs, Vsmax of adsorbates toward: O1, O1-O2, O2, ...
Forming catads:
Step1: Align X-D and A-X of cat and ads.
Step2: Using secondary axis of adsorbates to rotate horizontally to in_o3 vectors positions. 
Step3: Using primary axis of adsorbates to rotate vertically to apply tilt. 

Also, version 2 gives 0º < angle(Ti-X, normal o3) < 10º

Bind CO, Two sides: C and O;
    From previous works: its stable config is (adsorbate main axis) perp to(TiO3 plane)
    apply tilt phi=15°=pi/12 for a small perturbation. Check whether opt would push it back

Bind NH3, Two configs: eclipsed and staggered
    From previous works: its stable config is (adsorbate main axis) perp to(TiO3 plane)
    phi=15°=pi/12 for a small perturbation. Check whether opt would push it back

Bind H2O, One configs: 
    From previous works: its stable config is (adsorbate main axis)para to(TiO3 plane), face-lay
    But its Vsmin is at bottom. Using phi=60° for face-lay. 

Bind H2S, One config: 
    From pervious works: its stable config is (adsorbate main axis)para to(TiO3 plane), face-lay
    which is in consistent with its Vsmin

"""






"""
Store infos of adsorbate in a dict of lists:
[adsorbate, binding config, vsmin, vsmax, main axis, secondary axis, center atom index]
"""


class Ligand:
    def default_lig_config(self):
        """
        set Vsmin - Donor atom as the primary axis
        set Vsmax - Donor atom as the secondary axis
        """
        from lzn.atoms_construction.my_linear_algebra import normalize
        from lzn.atoms_construction.known_ligands import good_secondary_axis

        primary_axis = normalize(self.lig_atoms.positions[self.vsmin_at_idx] - self.vsmin_coord)
        secondary_axis = normalize(self.vsmax_coord - self.lig_atoms.positions[self.vsmin_at_idx])
        secondary_axis = good_secondary_axis(primary_axis, secondary_axis)
        self.lig_dict = {
            "c": [
                self.ligand_name,
                "c",
                0.,
                self.lig_vsmin[2],
                self.lig_vsmax[2],
                primary_axis,
                secondary_axis,
                self.vsmin_at_idx,
            ]
        }
        return [primary_axis, secondary_axis, self.vsmin_at_idx]

    def __init__(self, ligand_name, lig_hs95path = None):
        self.ligand_name = ligand_name.lower()
        if self.ligand_name in known_ligands.dict_lig.keys():
            self.lig_dict, self.lig_atoms = known_ligands.dict_lig[self.ligand_name]()
        else:
            if lig_hs95path is None:
                raise TypeError(f"{ligand_name} is not a predefined Ligand. lig_hs95path must be provided")
            self.lig_atoms = parse_standard_orientation_to_atoms(lig_hs95path)
            self.lig_vsmin = parse_vs_tolist(lig_hs95path, vs="min")[0]
            self.lig_vsmax = parse_vs_tolist(lig_hs95path, vs="max")[0]
            self.vsmin_at_idx = self.lig_vsmin[-1]
            self.vsmin_coord = self.lig_vsmin[1]
            self.vsmax_coord = self.lig_vsmax[1]
            self.default_lig_config()

    def simple_config(self):
        """
        Use one configuration for the ligands. Can be applied to CO and NH3 for simplification
        """
        for key_ in self.lig_dict.keys():
            if len(self.lig_dict[key_]) == 8:
                self.lig_dict[key_].pop(-2)
                self.lig_dict[key_][2] = 0.
            else:
                print("Already simple")

    # def apply_rot(self, rot = 2):
    #     """
    #     apply n-fold rotation to the simple configured ligand
    #         - rot: default C2
    #     """
    #     if len(self.lig_dict[key_]) == 8:



class Adsorption(Ligand):
    from pathlib import Path
    def __init__(self, ligand_name, tio2_hs95path, save_to_dir = Path.cwd(), **lig_hs95path):
        super().__init__(ligand_name, **lig_hs95path)
        self.tio2_hs95path = to_path(tio2_hs95path)
        self.save_to_dir = to_path(save_to_dir)

    def default_lig_config(self):
        super().default_lig_config()

    def simple_config(self, rot = "C1"):
        super().simple_config()

    # General binding algorithm for all adsorbates on all sites
    def build(self, align_to="tix"):
        """
        Read infos of cluster, using predefined functions.
        And accept infos of adsorbate by a list, which is written by later defined functions.
        Generate binding structures, and save to xyz file.
        Write Vs's in xyz comment line. Use Vsmin of O in cluster which adsorbate tilt of face toward.
        Kid of function calc_o3_plane_directions.
        l_adsorbate = [adsorbate, binding config name, tilt angle, vsmin, vsmax, primary axis, secondary axis, center atom index]
            - self.tio2_hs95path: cluster .out file
            - l_adsorbate: list of adsorbate info, receive from function bind_xxx
            - ads_atoms: list of ads info, receive from function bind_xxx
            - save_dir: directory to save the .xyz file
            - align_to: "tix" or "o3 normal"
        """
        import numpy as np
        from ase.io import write
        from lzn.atoms_construction.my_linear_algebra import robust_align, normalize
        from scipy.spatial.transform import Rotation

        more_config = True
        n_keys = 0
        for key_ in self.lig_dict:
            l_adsorbate = self.lig_dict[key_]
            n_keys += 1

            # Load adsorbate, use readable names as var
            if len(l_adsorbate) == 7:
                print("Use only one axis => One config for each site. ")
                more_config = False

            elif len(l_adsorbate) != 8:
                raise ValueError(
                    "len(l_adsorbate) must be 7/8, which are corresponding to one/two axis (one/more config for each site)")

            ads_titl_angle = l_adsorbate[2]  # Usage: calc binding
            ads_vsmin = l_adsorbate[3]  # Usage: write xyz comment
            ads_vsmax = l_adsorbate[4]  # Usage: write xyz comment

            ads_prim = l_adsorbate[5]  # Usage: calc binding
            if more_config:
                ads_secd = l_adsorbate[6]  # Usage: calc binding
            ads_center_pos = self.lig_atoms.positions[l_adsorbate[-1]]  # Center atom position (the atom which bind to Ti)

            # Read cluster
            cat_dict, cat_atoms = calc_o3_plane_directions(self.tio2_hs95path)
            cat_ino3_ll = cat_dict["dir_ino3"]  # list of lists: 6 vectors in O3 plane for each sites

            cat_vsmin_ll = cat_dict["cat_vsmin"]  # Usage: write xyz comment

            # filename infos
            fn_cat = self.tio2_hs95path.name
            if fn_cat[-8:] != "S_hs.out":
                raise SystemExit("Tore MSPA: Input file suffix should be S_hs.out => QUIT")
            fn_cat = fn_cat[:-8]
            fn_dict = parse_lzn_stem(fn_cat)

            fn_dict["ads"] = l_adsorbate[0]
            fn_dict["conf"] = l_adsorbate[1]

            ### Each binding site
            for site_, out_ in enumerate(cat_dict["out"]):  # Site number and Site infos
                catx_pos = out_[1]  # Catalyst Vsmax position
                catx_index = out_[0]  # Catalyst Vsmax index
                cat_vsmax = out_[2]  # Catalyst Vsmax value

                catti_pos = cat_atoms.positions[out_[-1]]
                match align_to:  # Choose catalyst's primary axis
                    case "tix":
                        cat_prim_vec = normalize(catx_pos - catti_pos)  # Catalyst primary vector
                    case "o3 normal":
                        o01_vec = cat_atoms.positions[out_[-1][1] - out_[-1][0]]
                        o02_vec = cat_atoms.positions[out_[-1][2] - out_[-1][0]]
                        cat_prim_vec = normalize(np.cross(o01_vec, o02_vec))
                    case _:
                        raise ValueError(f"Invalid align_to argument {align_to}")

                cat_sec_vec_l = cat_ino3_ll[site_]  # Catalyst secondary vector

                ### One config case:
                if not more_config:
                    cat_main_vec = cat_prim_vec
                    # Align adsorbate to the main axis, and center to 0.
                    ads_atoms_ = self.lig_atoms.copy()
                    rot_matrix = robust_align(ads_prim, cat_main_vec)
                    ads_atoms_.positions = (ads_atoms_.positions - ads_center_pos) @ rot_matrix.T
                    ## Translation
                    t_vec = catx_pos  # since ads already centered in previous rotation step, translation vector = position of X on cat
                    ads_atoms_.positions = ads_atoms_.positions + t_vec
                    # Format Comment and Filename
                    comment = f"#Written by construct_catads_adv# {cat_vsmax}; {ads_vsmin}"
                    save_name = f"{fn_dict['structure_name']}_{fn_dict['ads']}_site{catx_index:02d}_{fn_dict['conf']}.xyz"

                    # Construct cat_ads, and save to xyz
                    atoms_tot = cat_atoms + ads_atoms_
                    save_path = self.save_to_dir / save_name
                    write(save_path, atoms_tot, comment=comment)
                    conf_ = 0  # For print statement at the end

                ### More config case:
                else:
                    ### Each binding config
                    for conf_, cat_sec_vec in enumerate(cat_sec_vec_l):  # Conf number and Conf vector
                        ## Rotation1: align adsorbate to cluster primary axis with tilt
                        # (Can be solved using vector calc, but the result is too messy. Just use scipy Rotation for clarity)
                        tilt_angle = ads_titl_angle
                        if not np.isclose(tilt_angle, 0.):
                            rot_axis = normalize(np.cross(cat_prim_vec,
                                                          cat_sec_vec))  # It's safe to use cross product directly, since Ti-X is roughly perpendicular to O3-plane
                            rot_matrix = Rotation.from_rotvec(tilt_angle * rot_axis).as_matrix()
                            cat_main_vec = cat_prim_vec @ rot_matrix.T  # This vector will be the axis adsorbate to align with
                        else:  # scipy Rotation is not stable for very little tilt: (tilt_angle * rot_axis)->0
                            cat_main_vec = cat_prim_vec

                        # Align adsorbate to the main axis, and center to 0. Don't forget to rotate axes as well
                        ads_atoms_ = self.lig_atoms.copy()
                        rot_matrix = robust_align(ads_prim, cat_main_vec)
                        ads_atoms_.positions = (ads_atoms_.positions - ads_center_pos) @ rot_matrix.T
                        ads_secd_ = ads_secd @ rot_matrix.T

                        ## Rotation2: Move adsorbate face toward in_O3 vec directions
                        # Rotate around primary axis (absorbate primary axis has been already aligned with cluster's main axis)
                        # To move the secondary axis toward catalyst secondary axis
                        # Geometrical analysis -> align adsorbate secondary axis to vec_target = sec - dot(sec, prim)*prim/prim**2
                        target_vec = cat_sec_vec - np.dot(cat_sec_vec,
                                                          cat_main_vec) * cat_main_vec  # all the vectors are noramlized (CHECK!!!)
                        rot_matrix = robust_align(ads_secd_, target_vec)
                        ads_atoms_.positions = ads_atoms_.positions @ rot_matrix.T  # adsorbate's center atom is already centered

                        ## Translation
                        t_vec = catx_pos  # since ads already centered in previous rotation step, translation vector = position of X on cat
                        ads_atoms_.positions = ads_atoms_.positions + t_vec

                        # Format Comment and Filename
                        cat_vsmin = cat_vsmin_ll[site_][conf_]
                        comment = (f"#Written by construct_catads_adv# "
                                   f"{cat_vsmax} {cat_vsmin[0]} {cat_vsmin[1]}; "
                                   f"{ads_vsmin} {ads_vsmax}")

                        save_name = (f"{fn_dict['structure_name']}_{fn_dict['ads']}_"
                                     f"site{catx_index:02d}_{fn_dict['conf']}{conf_}.xyz")

                        # Construct cat_ads, and save to xyz
                        atoms_tot = cat_atoms + ads_atoms_
                        save_path = self.save_to_dir / save_name
                        write(save_path, atoms_tot, comment=comment)

        print(f"Written {site_ + 1} sites, with {n_keys} * {conf_ + 1} configurations for each, {n_keys * (site_ + 1) * (conf_ + 1)} in total.")
        return 0



if __name__ == "__main__":
    cluster = "/home/zining/TiO2_organized/05_mspa_clusters/tio12CDa_opt-pbe0_def2tzvpS_hs.out"
    save_to_dir = "~/tmp"

