# Wood Stiffness Homogenization

This code computes the multiscale stiffness of native, polymer-infiltrated, and delignified + infiltrated wood. It works from the cell-wall constituents up to clearwood. It is a procedural port of `02_stiffness/hom_TranspWood_model.m`, with a lumped cell wall: the middle lamella (ML) and S2 layer are not separated.

RVE chain: `pn → cel → cw → EWuc/LWuc → vesselwood → ringwood → ray → clearwood`

## References

[1] M. Königsberger, S. Scheiner, O. Lahayne, J. Schindler, D. Nuvoli, A. Mariani, J. Füssl:
*The micromechanics of transparent wood: Insights from multiscale modeling.*
Submitted to Composites Part B, 2026. (Model implemented here.)

Affiliations:
- Institute for Mechanics of Materials and Structures, TU Wien, Karlsplatz 13/202, 1040 Vienna, Austria
- Department of Chemical, Physical, Mathematical, and Natural Sciences, University of Sassari, Via Vienna 2, 07100 Sassari, Italy

[2] J. Schindler, M. Königsberger, L. Zelaya-Lainez, A. Mariani, D. Nuvoli, L. A. Berglund, J. Füssl:
*Nanoindentation of Wood Cell Walls: Effects of Delignification and Polymer Infiltration.*
Composites Science and Technology. (Cell-wall swelling; see in particular the Supplementary Information.)

## Files

| file | content |
|---|---|
| `homogenization_tw.py` | the whole model (tensors, Hill tensors, Mori–Tanaka / self-consistent schemes, volume fractions, RVE chain), plus input-file and batch handling and the command-line entry |
| `gui.py` | tkinter GUI |
| `woods.json` | species data: MFA, cell fractions, chemistry, density, EW/LW geometry and crystallinity CI, with references (converted from `01_data/data_Wood.xlsx`) |
| `polymers.json` | infiltration polymers (E, ν) |
| `phases.json` | cell-wall constituents: crystalline cellulose variants, amorphous cellulose, hemicellulose, lignin, densities |
| `examples/birch_batch.json` | example input file |
| `test_vs_matlab.py` | compares every RVE level with the MATLAB reference in `bench_matlab.mat` |
| `matlab/export_benchmark.m` | regenerates `bench_matlab.mat` from the MATLAB model; this is the only file that points outside the folder, to `02_stiffness` |

## Install

```bash
cd <PATH to folder with pyproject.toml>
pip install .
```

## Usage

### CLI

The installation will make available a `aitw-wood-stiffness-homogenization` command line interface.

- Run `aitw-wood-stiffness-homogenization --help` to see the available commands.
- Run `aitw-wood-stiffness-homogenization run --help` to see all available options.
- Run `aitw-wood-stiffness-homogenization run JSON_FILE` to run an homogenization calculation.

Example for birch microstructure generation:

```bash
aitw-wood-stiffness-homogenization run examples/birch_batch.json
```

#### Tab autocompletion

Enabling tab autocompletion https://click.palletsprojects.com/en/stable/shell-completion/

E.G for `bash` run the command

```bash
eval "$(_AITW_WOOD_STIFFNESS_HOMOGENIZATION_COMPLETE=bash_source aitw-wood-stiffness-homogenization)"
```

You can also add it to either `~/.bashrc` or, if you are using a virtual environment, to `bin/activate` of the virtual environment to avoid running the command for every new shell.

### GUI

**Start with**: `aitw-wood-stiffness-homogenization gui`

- Pick a wood. All fields are prefilled from `woods.json` and stay editable. Density and moisture are always shown; everything else is under *Advanced*.
- Pick a polymer, or type your own E and ν, and set the interface compliance and swelling.
- Every numeric field takes one value, a list (`550, 650`) or a range (`500:50:800`). All combinations are computed.
- Native, infiltrated and delignified + infiltrated are always computed. They are plotted in brown, blue and green.
- The plot depends on the batch:
  - nothing varied: bar chart of the moduli;
  - one input varied: x-y plot;
  - several inputs varied: choose the x-axis, and pick the values of the other varied inputs from dropdowns.
- Click a table row to see its volume fractions and 6×6 stiffness tensors.
- *Save input…* writes the current fields as an input file, and *Load input…* reads one back. *Save CSV* writes the results.

### Programmatically

```python
import json
from wood_stiffness_homogenization import homogeneization_tw as h

b = h.run_batch(json.load(open('examples/birch_batch.json')))   # or a dict built in code
h.write_csv(b, 'results.csv')
for d, r, i in b['rows']:              # one row per state, combination and MFA
    print(h.STATES[d['state']], d['density'], d['moisture'], d['E_poly'], r['E_L'][i])

p = h.defaults('Birch')                # single run with the full parameter dict
p.update(state=2, polymer='HEMA', MFA=[0, 10, 20])
r = h.run(p)                           # r['C_clearwood'], r['E_L'], r['levels'][i]['cw'], r['vol'], ...
```

## Input file

An input file is a JSON file. It uses the same keys as the GUI, so a file saved from the GUI can be run on the command line and the other way round. **Missing keys take the species defaults** from `woods.json`, so a minimal file is just `{"wood": "Spruce"}`.

Every numeric key takes a number, a list `[a, b, c]` or a range string `"start:step:stop"` (stop included). All combinations of all multi-valued keys are computed. Native wood ignores the polymer keys, so it is computed only once for each combination of the other inputs.

Example ([examples/birch_batch.json](examples/birch_batch.json): 4 densities × 2 moistures × 3 polymer moduli, with native computed once per density and moisture; 56 runs in total):

```json
{
 "wood": "Birch",
 "states": [1, 2, 3],
 "polymer": "HEMA",
 "density": "550:50:700",
 "moisture": [0, 12],
 "E_poly": [2.2, 4.5, 7.5],
 "nu_poly": 0.4,
 "MFA": [10],
 "swelling": true,
 "swelling_pct": 10,
 "output": "birch_batch_results.csv"
}
```

| key | meaning | default |
|---|---|---|
| `wood` | species name in `woods.json` | Birch |
| `states` | 1 native, 2 polymer infiltrated, 3 delignified + infiltrated | [1, 2, 3] |
| `polymer` | name in `polymers.json` (sets `E_poly`, `nu_poly` unless these are given) | HEMA |
| `E_poly`, `nu_poly` | polymer Young's modulus [GPa] and Poisson's ratio | from `polymer` |
| `density` | native wood density [kg/m³] | species mean |
| `moisture` | moisture content [% dry mass] | 12 |
| `MFA` | microfibril angles [°]; a list evaluated within each run, not a batch dimension | species mean |
| `CI` | mass-based crystallinity of cellulose | species |
| `cellulose`, `hemicellulose`, `lignin`, `extractives` | dry-mass fractions of the cell wall | species |
| `vessel_pct`, `ray_pct` | vessel and ray volume fractions [%] | species |
| `L_EW`, `W2_EW`, `L_LW`, `W2_LW` | earlywood/latewood lumen width and double wall thickness [µm]; `null` for hardwoods | species |
| `IF_alpha`, `IF_beta` | fibril–matrix interface compliance, tangential and normal [1/GPa] | 0, 0 |
| `swelling`, `swelling_pct` | wall swelling on/off, and increase of the wall thickness [%] | false, 10 |
| `cellulose_material` | `Dri2014_TI`, `IMWS` or `Dri` | Dri2014_TI |
| `cell_aspect_ratio`, `superellipse_n` | fibre-lumen slenderness, and lumen shape exponent (the ray-cell lumen is fixed at 1/5) | 0.01, 4 |
| `n_families`, `microfibril_slenderness`, `tolerance` | numerical settings | 20, 1e-20, 8 |
| `output` | result CSV, relative to the input file (command line only) | `<input>_results.csv` |

Keys starting with `_` are ignored, so you can use them for comments.

**Output CSV:** one line per state, combination and MFA. It contains all inputs, E_R/T/L, G_TL/RL/RT, the four Poisson's ratios, the cell-wall indentation modulus, and the upper triangles of the clearwood and cell-wall stiffness tensors (`Cclw_ij`, `Ccellwall_ij`, Kelvin–Mandel notation, GPa).

## Notes

- **Checked against MATLAB:** every RVE level agrees with the MATLAB model to about 1e-10 (`python test_vs_matlab.py`). The test covers Birch, Spruce and Pine; the native, infiltrated and delignified + infiltrated states; HEMA and PMMA; MFA sweeps; and interface compliance.
- **Ray orientation:** rays run radially, 90° from L (orientation `[0, pi/2]`; angles are in radians). Earlier versions of the MATLAB model used `[0, 90]`, which is 90 rad = 63.4° from L. That underestimated E_R (Birch with the former 1/100 ray lumen: 1.94 instead of 3.29 GPa) and produced a spurious C15 coupling. Fixed in Python and in `hom_TranspWood_model.m` / `hom_TranspWood_final.m`.
- **Non-symmetric Mori–Tanaka results:** with several orientation families the result is slightly non-symmetric. As in MATLAB, the raw tensor is reported, and its symmetric part is passed to the next scale.
- **Ray-cell lumen:** a spheroid with aspect ratio 1/5 along the ray axis. Fibre lumens use `cell_aspect_ratio` (1/100). The ray itself is infinitely long in the clearwood RVE. Compared with 1/100, this lowers E_R by about 10 % in native wood and 1.5–3 % in infiltrated wood; the other moduli change by less than 1.5 %.
- **Swelling** (ML and S2 lumped) follows the cell-wall swelling approach of Schindler et al. [2], see in particular their Supplementary Information, as implemented in `comp_volfrac_TranspWood`: the cell-centre distance stays at its native value, so the wall grows into the lumen.
  - Square fibre and ray cells: wall/half-width = 1-sqrt(porosity), scaled by (1+s).
  - Softwood EW/LW cells: the double wall 2W becomes 2W(1+s) at fixed outer cell size. The EW share stays fixed.
  - Vessels are unchanged.
  - The solid constituents keep their volume. The extra wall volume is added to the wall pore phase, which is polymer in states 2 and 3 and empty in state 1.
  - The balance is exact: pore' = (pore + k - 1)/k and x' = x/k, where k is the ratio of swollen to native wall volume. MATLAB instead adds `added/f_new` and renormalises, which is only correct to first order.
- **Not ported:** the `phase_control` study options of the MATLAB model.
