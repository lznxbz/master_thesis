from pathlib import Path

hartree_to_eV = 27.211399
hartree_to_kJdmol = 2625.5002
eV_to_hartree = 0.036749304951
eV_to_kJdmol = 96.485307499



def to_path(arg) -> Path:
    """
    Convert str to Path, accept existing Path, else raise error.
        - arg: str or Path, the file path to convert
    """
    if isinstance(arg, str):
        return Path(arg).expanduser().resolve()
    elif isinstance(arg, Path):
        return arg.expanduser().resolve()
    else:
        raise TypeError(f"Expected str or Path, got {type(arg).__name__}")



def batch_do(
        fxn: callable,
        in_dir: Path | str,
        matching: str | list[str],
        out_dir: Path | str | bool = False,
        **kwargs,
    ) -> list:
    """
    in_dir + matching -> in_path. Then deliver in_path, out_dir and **kwargs into fxn.
    And return list of fxn returned values.
    Parse in_dir and matching, deliver specific file path into fxn;Also accept **kwargs.
    But it doesn't generate output filename, nor save output files.
        - fxn: callable[Path, Path, dict]
        - in_dir: Path or str
        - matching: str or list[str], matching case for input files
        - out_dir: Path or str
        - **kwargs
    """
    in_dir = to_path(in_dir)
    if out_dir != False:
        if out_dir == "":
            out_dir = in_dir
        else:
            out_dir = to_path(out_dir)
            if not out_dir.exists():
                create_dir = input(f"{out_dir} does not exist, create it? [y/n] ")
                if create_dir.lower() == "y":
                    out_dir.mkdir(parents=True, exist_ok=True)
                else:
                    print("No folder created. --> Exit.")
                    return 0

    # Find input files
    if isinstance(matching, str):
        in_list = sorted(in_dir.glob(matching))
    elif isinstance(matching, list):
        in_list = [in_dir / m for m in matching]
    else:
        raise TypeError(f"Expected str or list, got {type(matching).__name__}")

    print(f"{len(in_list)} matching files")
    out_list = []
    for f in in_list:
        if out_dir != False:
            out_results = fxn(f, out_dir, **kwargs)
        else:
            out_results = fxn(f, **kwargs)

        if out_results != None:
            out_list.append(out_results)

    return out_list






def parse_lzn_stem(s: str) -> dict:
    """
    Return dict with keywords: structure_name, opt, ads, cm, sp
    Parse the stem of a LZN xyz file to extract structure name, opt method, sp method, comments, and adsorbate information.
    The expected format of the stem is: {structure_name}_opt-{opt_method}_sp-{sp_method}_cm-{comments}_ads-{adsorbate}
    where the parts after structure_name are optional and can be in any order.
    If a part is missing, the corresponding value in the returned dict will be None. 
        - s: the input string to parse
    """
    import re

    if isinstance(s, Path):
        s = s.stem

    keywords = ["opt", "ads", "cm", "sp"]
    pattern = r"_(opt|ads|cm|sp)-"

    # find all matches of the pattern
    matches = list(re.finditer(pattern, s))

    if not matches:
        info_dict = {
            "structure_name": s,
        }
    else:
        info_dict = {}
        # The first part before the first keyword is the structure name
        first_match = matches[0]
        info_dict["structure_name"] = s[:first_match.start()].rstrip("_-")
        # Extract values for each keyword
        for i, match in enumerate(matches):
            key = match.group(1) # the keyword in the parentheses r"_(opt|ads|cm|sp)-"

            start = match.end()
            if i + 1 < len(matches):
                end = matches[i + 1].start()
            else:
                end = len(s)

            value = s[start:end].rstrip("_-")
            info_dict[key] = value if value else None

    return info_dict



def write_lzn_stem(dict_: dict) -> str:
    """
    The inverse of parse_lzn_stem. Receive the information dictionary, write to filename stem.
    """
    stem = dict_["structure_name"]
    for key in dict_.keys():
        if key == "structure_name":
            continue
        stem += f"_{key}-{dict_[key]}"

    return stem





def check_mspa(s: str) -> tuple[bool, str]:
    """
    Check if the structure is for MSPA (ends with "S_hsM"), which record the X sites
    and return the cleaned string without "S_hsM".
    Particularly, used for the step of binding sub and ads to subads.
    subS_hsM + adsS_hsM -> sub_ads-ads
        - s: the input string to check
    """
    if s[-5:] == "S_hsM":
        mspa = True
        s_clean = s[:-5]
    else:
        mspa = False
        s_clean = s
    return mspa, s_clean




class LZNxyz:
    def __init__(
            self, 
            structure_name: str, 
            opt_method: str,
            file_path: Path,
            mspa: bool = False,
            energy_hartree: float = 0.0, 
            from_structure: str = None, 
            sp_method: str = None, 
            comments: str = None, 
            adsorbate: str = None,
            **kwargs
        ):
        if file_path.suffix != ".xyz":
            raise ValueError(f"File {file_path} is not an xyz file")
        self.structure_name = structure_name
        self.file_path = file_path
        self.opt_method = opt_method
        self.mspa = mspa
        self.sp_method = sp_method
        self.comments = comments
        self.adsorbate = adsorbate
        self.energy_hartree = energy_hartree
        self.from_structure = from_structure

    def lzn_xyz_stem(self) -> str:
        # cond = (self.from_structure == None) or (self.from_structure == self.structure_name)
        # if cond:
        #     return f"{self.structure_name}_opt-{self.opt_method}"
        # else:
        #     return f"{self.structure_name}_opt-{self.opt_method}_cm-{self.from_structure}"
        return f"{self.structure_name}_opt-{self.opt_method}"

    def lzn_xyz_comment(self) -> str:
        comment = f"E={self.energy_hartree} hartree"
        if self.from_structure:
            comment = f"from {self.from_structure}; {comment}"
        return comment

    def energy_in_eV(self) -> float:
        return self.energy_hartree * hartree_to_eV




def read_lzn_xyz(file_path: Path | str, read_xyz_energy: bool = False) -> LZNxyz:
    """
    Read a LZN xyz file and return a LZNxyz object containing the structure and energy information.
    Parse the filename to extract structure name, opt method, sp method, comments, and adsorbate information.
    Also check if the structure is for MSPA based on the filename.
        - file_path: the path to the xyz file to read
    """
    import re
    path = to_path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"File {file_path} does not exist")

    # Parse filename
    s = path.stem
    mspa, s_clean = check_mspa(s)
    info_dict = parse_lzn_stem(s_clean)

    if mspa or (not read_xyz_energy):
        energy_hartree = 0.0
    else:
        # Read energy from comments line
        with path.open() as f:
            lines = f.readlines()
            if len(lines) < 2:
                raise ValueError("LZN xyz file must have at least 2 lines")
            comment_line = lines[1].strip()

        energy_match = re.search(r'E\s*=\s*([+-]?\d*\.?\d*)\s*([a-zA-Z]+)', comment_line)
        if not energy_match:
            print(f"Could not parse energy from line: {comment_line}")
            energy_hartree = 0.0
        else:
            energy_unit = energy_match.group(2)
            if energy_unit == "eV":
                energy_hartree = float(energy_match.group(1)) / hartree_to_eV
            elif energy_unit in ["hartree", "Hartree"]:
                energy_hartree = float(energy_match.group(1))
            else:
                print(f"Unrecognized energy unit: {energy_unit}")

    return LZNxyz(
        structure_name=info_dict["structure_name"], 
        opt_method=info_dict["opt"], 
        file_path=path, 
        mspa=mspa,
        energy_hartree=energy_hartree,
        sp_method=info_dict["sp"],
        comments=info_dict["cm"],
        adsorbate=info_dict["ads"]
    )






if __name__ == "__main__":
    test_lznxyz = read_lzn_xyz("~/TiO2_organized/06_cluster_ads_xyz/tio12CDa_opt-omol25_low_ads-nh3_site36.xyz")
    print(test_lznxyz.structure_name)
    print(test_lznxyz.opt_method)
    print(test_lznxyz.sp_method)
    print(test_lznxyz.comments)
    print(test_lznxyz.adsorbate)




