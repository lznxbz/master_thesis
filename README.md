## Repository Structure
```text
.
├── atoms_construction
│   ├── construct_complex_v3.py
│   ├── data_hs95_out
│   │   ├── coS_hs.out
│   │   ├── h2oS_hs.out
│   │   ├── h2sS_hs.out
│   │   ├── nh3S_hs.out
│   │   └── oh_radicalS_hs.out
│   ├── __init__.py
│   ├── known_ligands.py
│   ├── my_linear_algebra.py
│   ├── parse_hs95_out.py
│   └── summary_all_possible_config.py
├── data_analysis
│   ├── check_g16_output.py
│   ├── extract_energies.py
│   ├── find_symmetry_equivalent_sites.py
│   ├── __init__.py
│   └── rmsd_min.py
├── file_convert
│   ├── aseopt_traj_to_xyz.py
│   ├── chk_to_g16orb_gjf.py
│   ├── g16opt_log_to_xyz.py
│   ├── __init__.py
│   └── write_gjf.py
├── general.py
├── __init__.py
├── misc
│   ├── check_consistency_dft_mlip.py
│   └── check_consistency_with_symmetries.py
└── run_calc
    ├── calc_eos.py
    ├── __init__.py
    └── run_mlip_v1.py

```
## 1. Construct Adsorption Complex Configurations
This module was developed specifically for $TiO_2$ clusters and has only been tested on the **12a** and **5a** $TiO_2$ clusters from *Nanoscale, 2017, 9, 7143*.
### File Specifications
 * **Input files:** HS95 output files (xxxS_hs.out) for the cluster and ligands. (Predefined ligands include: h2o, h2s, co, nh3, and oh_radical).
 * **Output files:** Structural .xyz files for the adsorption complexes (total files generated = number of sites $\times$ number of configurations).
### Usage and APIs
```python
# Import the Adsorption class module
from lzn.atoms_construction.construct_complex_v3 import Adsorption

# Initialize the class by specifying the ligand name (chemical formula string) and the path to the cluster hs95.out file (string or Path)
adsorption = Adsorption(ligand_name, tio2_path) 

# If the ligand is not predefined (i.e., not h2o, h2s, nh3, co, or oh_radical), 
# provide the path to the ligand using the 'lig_hs95path' keyword argument:
adsorption = Adsorption(ligand_name, tio2_path, lig_hs95path=ligand_path) 

# The default save directory is the current working directory (CWD).
# To specify a custom output directory, use the 'save_to_dir' keyword argument:
adsorption = Adsorption(ligand_name, tio2_path, save_to_dir=out_dir, lig_hs95path=ligand_path) 

# Alternatively, the output directory can be assigned via the object's property:
adsorption.save_to_dir = to_path(save_to_dir)

# To generate only one configuration per site, use the simple_config() method:
adsorption.simple_config()

# To apply the default universal binding configurations for predefined ligands, use the default_lig_config() method:
adsorption.default_lig_config()

# For a more customized binding configuration, modify the 'lig_dict' property.
# Format: {"config name": [ligand_name, "config name", tilt_angle, lig_vsmin value, lig_vsmax value, primary_axis, secondary_axis, vsmin_at_atom index]}
adsorption.lig_dict = xxxxx

# Once setup is complete, build the adsorption complex structures using the build() method.
# The 'align_to' argument accepts two options: "tix" (default) or "o3 normal"
# This will generate (number of sites * number of configurations) of structural .xyz files under adsorption.save_to_dir
adsorption.build() 

```
### Example
```python
from lzn.atoms_construction.construct_complex_v3 import Adsorption
from lzn.general import to_path

tio2_path = to_path('/path/to/tio2 cluster/hs95/.out file')
adsorption = Adsorption(ligand_name, tio2_path) 
adsorption.build() 

```
## 2. Run MLIP Optimization
### File Specifications
 * **Input files:** .xyz structural files.
 * **Output files:** ASE binary trajectory files (.traj) and energy log files (.log, in eV) tracking each step.
### Usage
```python
from lzn.run_calc.run_mlip_v1 import *

in_dir = "/dir/to/input/xyz/files/"  # Path or string
matching = "matching filenames"       # Supports wildcard (*) and regular expressions
out_dir = "/dir/to/save/output/files/" # Path or string

# Load the calculator and apply a shift if the modal is "omol25"
calc = shift_calc(SevenNetCalculator(model="7net-omni", modal="matpes_r2scan")) 

batch_do(
    run_mlip_opt,
    in_dir,
    matching,
    out_dir,
    calc=calc,
    out_xyz_dir="/dir/to/save/xyz/",  # Optional
    fmax=0.01,                        # Optional
) 

```
### Example
```python
from lzn.run_calc.run_mlip_v1 import *
from lzn.general import to_path

in_dir = "/home3/zining/backup/organized/01_cluster_article_xyz/"
matching = "*"
out_dir = to_path(__file__).parent

# Load the calculator and apply a shift if the modal is "omol25"
calc = shift_calc(SevenNetCalculator(model="7net-omni", modal="matpes_r2scan")) 

batch_do(
    run_mlip_opt,
    in_dir,
    matching,
    out_dir,
    calc=calc,
)

```
## 3. Format Converter
Export optimized structures to .xyz format using the modules below. These can be used in combination with lzn.general.batch_do.
 * lzn.file_convert.aseopt_traj_to_xyz.asetraj_to_xyz
 * lzn.file_convert.g16opt_log_to_xyz.log_to_xyz
### Example
```python
batch_do(
    log_to_xyz, 
    '/input/dir/',
    'matching',
    '/output/dir/',
)

```
## 4. Extract Information
Extract calculations and structural data using batch_do(fxn, in_dir, matching) combined with the following functions:
 * **Check Gaussian 16 Output / Frequencies:** lzn.data_analysis.check_g16_output.check_imaginary_frequencies
 * **Extract ASE Energies:** lzn.data_analysis.extract_energies.extract_ase_energy
 * **Extract Gaussian 16 Energies:** lzn.data_analysis.extract_energies.extract_g16_energy
## 5. Calculate Minimum RMSD
Calculate the minimum Root-Mean-Square Deviation (RMSD) between atomic configurations:
```python
lzn.data_analysis.rmsd_min

```
### API Methods
```python
# Find minimum RMSD between two full atomic structures
find_rmsd_min(
    a: Atoms,
    b: Atoms,
)

# Find minimum RMSD for partial atomic structures based on character lists
find_rmsd_min_partial(
    a,
    b,
    char_l,
)

```

## More Examples
see /home3/zining/pylzn/test/