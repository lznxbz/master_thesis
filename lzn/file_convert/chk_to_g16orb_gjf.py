from lzn.general import *

def sp_chk_to_gjf(file_in: Path | str, output_dir: Path | str = None) -> str:
    """
    Write a gjf file for calculating MOs from a sp chk file.
        - file_in: input sp chk file
        - output_dir: directory to save gjf file, if None, do not write gjf file and only return the gjf head as string
    """
    # Generate gjf head
    gjf_head_list = ["%mem=12gb", "pop=min guess=(read,only) geom=allcheck"]

    file_in = to_path(file_in)
    if not file_in.is_file():
        raise ValueError(f"Input file {file_in} does not exist or is not a file.")
    chk_name = file_in.stem
    gjf_head_list.insert(1, f"%oldchk={file_in.resolve()}")

    info_dict = parse_lzn_stem(chk_name)
    structure_name = info_dict["structure_name"]
    opt_method = info_dict["opt"]
    sp_method = info_dict["sp"]

    if sp_method == "pbe0_def2tzvp":
        gjf_head_list.insert(2, "#pbe1pbe/def2tzvp iop(6/7=3) gfinput 6d 10f")
    else:
        raise ValueError(f"Unsupported sp method: {sp_method}")
    
    gjf_head = "\n".join(gjf_head_list) + "\n\n"

    # Write gjf file
    if output_dir is None:
        print("No ooutput.")
    else:
        output_dir = to_path(output_dir)
        gjf_name = f"{structure_name}_opt-{opt_method}S.gjf"
        gjf_path = output_dir / gjf_name
        gjf_path.write_text(gjf_head, encoding="utf-8")
        print(f"GJF file written to {gjf_path}")

    return gjf_head





if __name__ == "__main__":
    file_in = "~/TiO2_organized/04_cluster_g16_sp/tio12CDa_sp-pbe0_def2tzvp_opt-matpes_r2scan.chk"
    output_dir = "~/TiO2_organized/05_mspa_substrates"
    print(sp_chk_to_gjf(file_in, output_dir))



