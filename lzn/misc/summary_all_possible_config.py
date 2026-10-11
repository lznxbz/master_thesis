from lzn.data_analysis.rmsd_min import find_rmsd_min, find_rmsd_min_partial
from lzn.general import to_path

from ase.io import read
import matplotlib.pyplot as plt
import pickle

"""
compare energy and rmsd of adv configs and ref
both dE, RMSD < tol => only one stable config.
"""

dict_conf = {
    "co": ["Cb", "Ob"],
    "h2o": ["c"],
    "h2s": ["c"],
    "nh3": ["1HP", "2HP"],
}

dict_natoms = {
    "co": 2,
    "h2o": 3,
    "h2s": 3,
    "nh3": 4,
}

dict_dir = {
    "matpes_r2scan": "771_r2scan_opt",
    "omol25_low": "770_omol25_opt",
    "pbe0_def2svp": "772_svp_opt",
}

with open(
    "/home/zining/organized/05_mspa_clusters/tio12CDa_opt-pbe0_def2tzvp_parse.pkl",
    "rb"
) as f:
    # (x index, x coord, x val, NN type=Ti, NN index, nearest 3O index) from parse_tore_mspa_out
    ll = pickle.load(f)
    # {site: [tio3 and ads indices]}
    tio2_info_dict = {}
    for l in ll:
        tio2_info_dict[l[0]] = [l[-2]] + [i for i in l[-1]]


def calc_dE_rmsd(abs: str, opt: str, rmsd_option = "total"):
    """
    - abs: absorbate molecule name
    - opt: optimization method
    - rmsd_option: rmsd calculation method, "total" OR "partial".
        "partial" only compare (TiO3 slab + adsorbate)
    Return: dict{site_: [(conf0_E, conf0_RMSD)], [()], ...},
            Emin_name["confx", ...]
    """
    ref_dir = to_path(f"~/organized/{dict_dir[opt]}")
    ref_filename = lambda site_: f"tio12CDa_{abs.lower()}_site{site_:02d}*_opt-{opt}.traj"

    all_dir = to_path("~/organized/07_adv_mlip_catads_opt")
    all_filename = lambda site_, conf_: f"tio12CDa_{abs.upper()}_site{site_:02d}_{conf_}_opt-{opt}.traj"

    dE_rmsd_dict = {}
    Emin_l = []
    Emin_name_l = []
    for site_ in range(0, 22, 2): # iteration over all sites
        if site_ == 16:
            continue
        ref_file = list(ref_dir.glob(ref_filename(site_)))
        if len(ref_file) == 1:
            ref_file = ref_file[0]
        else:
            raise IndexError(f"{len(ref_file)} matching files detected. Expect 1")
        ref_atoms = read(ref_file)
        ref_energy = ref_atoms.get_potential_energy() # Reference energy

        conf_l = dict_conf[abs.lower()]
        dE_rmsd_l = []
        Emin = 0.
        for i_ in conf_l: # 1 or 2 types of configurations (e.g. CO: C-bind, O-bind)
            for conf_ in range(6): # iteration over all configurations
                conf_ = i_ + str(conf_)
                all_file = all_dir / all_filename(site_, conf_)
                all_atoms = read(all_file)

                E_conf = all_atoms.get_potential_energy()
                if E_conf < Emin:
                    Emin = E_conf
                    Emin_name = conf_
                dE = E_conf - ref_energy # Energy relative to reference
                match rmsd_option.lower():
                    case "total":
                        x_, x__, rmsd = find_rmsd_min(ref_atoms, all_atoms)
                    case "partial":
                        char_l = tio2_info_dict[site_] + list(range(-dict_natoms[abs]))
                        x_, x__, rmsd = find_rmsd_min_partial(ref_atoms, all_atoms, char_l)
                    case _:
                        raise ValueError("Available options: total, partial")
                dE_rmsd_l.append((dE, rmsd))

        dE_rmsd_dict[f"site{site_}"] = dE_rmsd_l
        Emin_l.append(Emin)
        Emin_name_l.append(Emin_name)

    return dE_rmsd_dict, Emin_name_l



def plot_results(data, Emin_names):
    fig, axes = plt.subplots(4, 3, figsize=(15, 16))
    axes = axes.flatten()  # Flatten to iterate easily

    for i, (title, values) in enumerate(data.items()):
        ax1 = axes[i]
        # 2. Extract the data points
        indices = list(range(len(values)))
        y_left = [v[0] for v in values]   # Energy
        y_right = [v[1] for v in values]  # RMSD

        # 3. Plot Left Y-Axis (Energy)
        color_a = 'tab:blue'
        ax1.set_xlabel('Index')
        ax1.set_ylabel('dE/eV', color=color_a, fontsize=10)
        ax1.scatter(indices, y_left, color=color_a, marker='o', label='dE')
        ax1.tick_params(axis='y', labelcolor=color_a)

        # 4. Create and Plot Right Y-Axis (RMSD)
        ax2 = ax1.twinx()  # Create the twin axis
        color_b = 'tab:red'
        ax2.set_ylabel('rmsd/Å', color=color_b, fontsize=10)
        ax2.scatter(indices, y_right, color=color_b, marker='o', label='rmsd')
        ax2.tick_params(axis='y', labelcolor=color_b)

        title = f"{title}; min={Emin_names[i]}"
        ax1.set_title(title, fontweight='bold', pad=10)
        ax1.grid(True, alpha=0.3)

    if len(data) < len(axes):
        fig.delaxes(axes[-1])

    fig.tight_layout()
    plt.show()

    return 0



def main():
    data, Emin_names = calc_dE_rmsd("nh3", "omol25_low", "partial")
    plot_results(data, Emin_names)
    
    return 0


if __name__ == "__main__":
    main()
