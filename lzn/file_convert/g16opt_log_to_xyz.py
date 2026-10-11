from lzn.general import *
from ase.io.gaussian import read_gaussian_out

# extract dft opt from .log to .xyz

def log_to_xyz(in_path, out_path, index = -1):
    # read coordination
    in_path, out_path = map(to_path, [])


    fcontent_i = open(in_path, "r")
    atoms_i = read_gaussian_out(fcontent_i, index=index)

    # write xyz
    atoms_i.write(out_path, format="xyz")

    return 0





def batch_log_to_xyz(dft_opt_folder,
                     save_folder,
                     file_name: str = "*.log",
                     index = -1):
    for file_i in list(dft_opt_folder.glob(file_name)):
        # xyz file name
        fname_i = file_i.stem.split("_")
        fstructure_i = fname_i[0]
        fopt = fname_i[1]
        if (fopt == "pbe" or fopt == "pbe0"):
            fmethod_i = fopt + "_def2svp"
        elif (fopt == "pbetzvp" or fopt == "pbe0tzvp"):
            fmethod_i = fopt[:-4] + '_def2' + fopt[-4:]
        elif (fopt == "pbesvp" or fopt == "pbe0svp"):
            fmethod_i = fopt[:-3] + '_def2' + fopt[-3:]
        else:
            fmethod_i = input(f"{file_i.stem}, No match method, Input: (functional_basis)")

        save_name = fstructure_i + "_" + fmethod_i + ".xyz"
        print(save_name)
        save_path = save_folder / save_name

        log_to_xyz(file_i, save_path, index=index)

    print("Done")
    return 1




if __name__ == "__main__":
    fxn = lambda x: x+1
    ls = list(range(5))
    f = map(fxn, ls)

    print(f[1])