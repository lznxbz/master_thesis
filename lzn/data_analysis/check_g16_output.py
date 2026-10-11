from lzn.general import to_path, batch_do
import cclib
from colorama import Fore, Back, Style



def check_imaginary_frequencies(filepath, *args):
    """
    check imaginary frequencies
    """
    stem = filepath.stem
    data = cclib.io.ccread(filepath)
    try:
        frequencies = data.vibfreqs
    except AttributeError:
        print(Back.YELLOW + f"Skipping {stem}. No vibfreqs available")
        print(Style.RESET_ALL, end="")
        return 1

    lowest_frequency = frequencies[0]

    if lowest_frequency < 0.:
        print(Fore.RED + f"{stem} Warning: imaginary frequencies {lowest_frequency} < 0")
        print(Style.RESET_ALL, end="")
    else:
        print(f"{stem}: {lowest_frequency}")

    return 0



if __name__ == "__main__":
    batch_do(
        check_imaginary_frequencies,
        "/home/zining/TiO2_organized/772_g16_opt",
        "tio12CDa_*_opt-pbe0_def2svp.log",
        "",
    )


