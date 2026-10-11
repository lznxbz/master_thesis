from lzn.atoms_construction.parse_hs95_out import *
from pathlib import Path

dir_hs95out = Path(__file__).parent / "data_hs95_out"


# Read adsorbates MSPA .out file
def get_ads_info(ads):
    ads_path = dir_hs95out / f"{ads.lower()}S_hs.out"
    # Read adsorbate. Two scenarios staggered and eclipsed
    ads_atoms = parse_standard_orientation_to_atoms(ads_path)
    ads_vsmin = parse_vs_tolist(ads_path, vs="min")
    ads_vsmax = parse_vs_tolist(ads_path)
    return ads_atoms, ads_vsmin, ads_vsmax

def good_secondary_axis(primary_axis, secondary_axis):
    """
    In following bind_ads functions, the input secondary axis is not necessarily perpendicular to primary axis.
    It is only a representative axis that reflect molecular dipole moment properties.
    This function get the good secondary axis which is the projection in the plane with primary axis as normal axis.
    """
    costheta = np.dot(normalize(secondary_axis), normalize(primary_axis))
    good_sec = secondary_axis
    if abs(costheta) > 0.995 :
        raise ZeroDivisionError("Primary axis and secondary axis are overlapped")
    elif abs(costheta) > 1E-3:
        good_sec = secondary_axis - costheta * primary_axis
    return good_sec



def bind_co():
    """
    Bind CO, Two sides: C and O;
    tilt=18°=pi/10 wrt Ti-X axis for a small perturbation. Check whether opt would push it back
    """
    ads_atoms, ads_vsmin, ads_vsmax = get_ads_info("co")

    # Using C-O vectors as primary axis. It is more precise than X-C vector,
    # since FEM used to find X -> precision loss.
    primary_axis1 = normalize(ads_atoms.positions[1] - ads_atoms.positions[0])
    primary_axis2 = - primary_axis1

    secondary_axis = np.array([1., 0., 0.]) # CO is linear in z-axis
    secondary_axis = good_secondary_axis(primary_axis1, secondary_axis)

    co_dict = {
        "Cb": ["CO", "Cb", np.pi/10, ads_vsmin[0][2], ads_vsmax[0][2], primary_axis1, secondary_axis, 0],
        "Ob": ["CO", "Ob", np.pi/10, ads_vsmin[1][2], ads_vsmax[0][2], primary_axis2, secondary_axis, 1],
    }
    return co_dict, ads_atoms

def bind_nh3():
    """
    Bind NH3, Two configs: eclipsed and staggered. <=> (1HP, 2HP) + (o1, o1o2, o2, ...)
    phi=18° wrt Ti-X axis for a small perturbation. Check whether opt will push it back.
    Denote: H3 = H3 plane = H3 center, depend on context
    """
    ads_atoms, ads_vsmin, ads_vsmax = get_ads_info("nh3")

    # Using N-H3 vector as primary axis. It is more precise than X-N vector, since FEM used to find X -> precision loss.
    primary_axis = normalize(np.average(ads_atoms.positions[1:], axis=0) - ads_atoms.positions[0])

    # Using one of H3-H vector as secondary axis
    secondary_axis1 = normalize(ads_atoms.positions[1] - np.average(ads_atoms.positions[1:], axis=0))
    secondary_axis1 = good_secondary_axis(primary_axis, secondary_axis1)
    secondary_axis2 = - secondary_axis1

    ads_dict = {
        "1HP": ["NH3", "1HP", np.pi/10, ads_vsmin[0][2], ads_vsmax[0][2], primary_axis, secondary_axis1, 0],
        "2HP": ["NH3", "2HP", np.pi/10, ads_vsmin[0][2], ads_vsmax[0][2], primary_axis, secondary_axis2, 0],
    }
    return ads_dict, ads_atoms

def bind_h2o():
    """
    Bind H2O, One config. phi=60º=pi/3.
    Denote: H2 = H2 line = H2 center, depend on context.
    """
    ads_atoms, ads_vsmin, ads_vsmax = get_ads_info("h2o")

    # Using O-H2 vector as primary axis. It is more precise than X-O vector, since FEM used to find X -> precision loss.
    primary_axis = normalize(np.average(ads_atoms.positions[1:], axis=0) - ads_atoms.positions[0])

    # Secondary axis perpendicular to H2O plane
    h2_vec = ads_atoms.positions[1] - ads_atoms.positions[2]
    secondary_axis = normalize(np.cross(h2_vec, primary_axis))
    # secondary_axis = normalize(np.average(ads_atoms.positions[1:], axis=0) - ads_atoms.positions[0])
    secondary_axis = good_secondary_axis(primary_axis, secondary_axis)

    ads_dict = {
        "only": ["H2O", "c", np.pi/3, ads_vsmin[0][2], ads_vsmax[0][2], primary_axis, secondary_axis, 0],
    }
    return ads_dict, ads_atoms

def bind_h2s():
    """
    Bind H2S, One config. phi=0.
    Denote: H2 = H2 line = H2 center, depend on context.
    """
    ads_atoms, ads_vsmin, ads_vsmax = get_ads_info("h2s")

    # Primary axis = X-S vector
    primary_axis = normalize(ads_atoms.positions[0] - ads_vsmin[0][1])

    # S-H2 vector can be used as Secondary axis, since large angle(X-S-H2)~100º
    secondary_axis = normalize(np.average(ads_atoms.positions[1:], axis=0) - ads_atoms.positions[0])
    secondary_axis = good_secondary_axis(primary_axis, secondary_axis)

    ads_dict = {
        "only": ["H2S", "c", 0., ads_vsmin[0][2], ads_vsmax[0][2], primary_axis, secondary_axis, 0],
    }
    return ads_dict, ads_atoms

def bind_oh_radical():
    """
    Bind OH, One config. phi=0.
    """
    ads_atoms, ads_vsmin, ads_vsmax = get_ads_info("oh_radical")

    # Primary axis = X-O vector
    o_index = ads_atoms.get_chemical_symbols().index("O")
    primary_axis = normalize(ads_atoms.positions[o_index] - ads_vsmin[0][1])

    # O-H vector can be used as Secondary axis, since large angle(X-S-H2)~83º
    secondary_axis = normalize(np.average(ads_atoms.positions[1:], axis=0) - ads_atoms.positions[0])
    secondary_axis = good_secondary_axis(primary_axis, secondary_axis)

    ads_dict = {
        "only": ["OH_radical", "c", 0., ads_vsmin[0][2], ads_vsmax[0][2], primary_axis, secondary_axis, 0],
    }
    return ads_dict, ads_atoms


dict_lig = {
    "co": bind_co,
    "h2o": bind_h2o,
    "nh3": bind_nh3,
    "h2s": bind_h2s,
    "oh_radical": bind_oh_radical,
}
