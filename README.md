# mecp2_paper

Analysis code and figures for *Local circuit connectivity enhances visual decoding and
visuomotor behavior* (Privat, Hansen, Pietri et al.).

The code compares wild-type and *mecp2*-mutant zebrafish larvae:

- **Figure 1:** spontaneous activity in the optic tectum (two-photon calcium imaging): event
  frequency and duration, pairwise correlations against a shuffled null model.
- **Figure 2:** tectal responses to looming stimuli of increasing contrast, and decoding of
  contrast from population activity.
- **Figure 3:** tectal responses to small prey-like spots at four azimuths, and decoding of
  azimuth. Supplementary figure: example responses and decoding error.
- **Figure 4 (E) and Supplementary Figure 2:** J-turns towards virtual prey of increasing UV
  intensity in a closed-loop virtual reality assay, and sensitivity (EC50) from bootstrapped
  sigmoid fits.

Figures 1–3 are produced in MATLAB, Figure 4 in Python. Statistics treat the larva as the
independent unit: per-larva values are compared with Wilcoxon rank-sum tests, or with a
mixed-design ANOVA (genotype between larvae, stimulus within larvae, Greenhouse–Geisser
corrected) for responses across stimulus levels.

## Quick start

```bash
git clone https://github.com/ElTinmar/mecp2_paper.git
cd mecp2_paper
conda env create -f environment.yml
conda activate mecp2
python download_data.py            # downloads, verifies and extracts the data into data/
python run_all.py data/DATA        # figures and tables in results/
```

MATLAB R2024a with the toolboxes listed under [Requirements](#requirements) must be installed
and on your `PATH`.

## Running

From the conda environment, with MATLAB on your `PATH`:

```bash
python run_all.py data/DATA                        # all figures, in parallel (64GB+ RAM)
python run_all.py data/DATA --jobs 2               # at most two at a time (32GB+ RAM)
python run_all.py data/DATA --figures F2 F4 --out results_v2
```

Each figure runs in its own process. Outputs go to `results/` (or `--out`), and each figure's
console output to `results/logs/`. Lower `--jobs` if memory or MATLAB licenses are limited.
Use `--matlab /path/to/matlab` if MATLAB is not on your `PATH`.


## Outputs

| Figure | Figures | Tables | Other |
| --- | --- | --- | --- |
| 1 | `fig1.pdf`, `fig1.svg` | `table_F1.txt` | `fig1.mat` |
| 2 | `fig2.pdf`, `fig2.svg` | `table_F2.txt` | `fig2.mat`, `SupMovie_F2.avi` |
| 3 | `fig3.pdf`, `fig3.svg`, `figS2.pdf`, `figS2.svg` | `table_F3.txt` | `fig3.mat`, `SupMovie_F3.avi` |
| 4 | `fig4.pdf`, `fig4.svg`, `figS2_F4.pdf`, `figS2_F4.svg` | `table_F4.txt` | `fig4.npz` |

The tables contain every statistic reported in the paper (cell counts, test statistics,
p-values, effect sizes).

## Requirements

**MATLAB R2024a** (figure export does not work with R2025a+) with:

- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox
- Image Processing Toolbox
- Curve Fitting Toolbox
- Computer Vision Toolbox

**Python**, from conda-forge:

```bash
conda env create -f environment.yml
conda activate mecp2
```

The figures use Arial. If it is not installed, matplotlib silently falls back to another font
(check with `pdffonts fig4.pdf`).

### On Linux

The libstdc++ version used by MATLAB can conflict with the system libraries. To fix it:

```bash
sudo apt install matlab-support
```

## Data

### Download

```bash
python download_data.py                       # into data/
python download_data.py --dest /path/to/data  # elsewhere
```

The data are a single ZIP archive. The script checks the archive against the SHA-256
recorded in `download_data.py`, and every file inside it (ZIP CRC-32) before extracting. A
corrupt download is renamed `.corrupt` and the script stops; run it again to retry.
Interrupted downloads resume. Use `--keep-zip` to keep the archive after extraction.

### Layout

The archive extracts to `DATA/`, the folder passed to `run_all.py`:

```
DATA/
├── Spontaneous/                                    Figure 1 (9 WT, 10 MeCP2 larvae)
│   ├── WT/
│   │   └── <date>_nacre_5dpf_<n>_reg_smooth_CELL_DATA.mat     deltaFoF, raster, movements, orig_position
│   └── MeCP2/
│       └── <date>_mecp2_5dpf_<n>_reg_smooth_CELL_DATA.mat
├── Stimulation_Contrast/                           Figure 2 (7 WT, 10 MeCP2 larvae)
│   ├── WT/
│   │   ├── <date>_6dpf_nacre_<nnn>.mat                        stimulation: paramContrast, timeContrast
│   │   └── <date>_6dpf_nacre_<nnn>_reg_smooth_CELL_DATA.mat   deltaFoF, raster, movements, avg_clean, cell
│   └── MeCP2/
│       └── <date>_6dpf_mecp2_<nnn>{.mat, _reg_smooth_CELL_DATA.mat}
├── Stimulation_Azimuth/                            Figure 3 (7 WT, 9 MeCP2 larvae)
│   ├── WT/
│   │   ├── <date>_GCaMP5_6dpf_<nn>_1-95Hz.mat                 stimulation: randomOrder, timeStart, timeStop
│   │   ├── <date>_GCaMP5_6dpf_<nn>_1-95Hz_reg_ALL_CELLS.mat   avg, cells
│   │   ├── <date>_GCaMP5_6dpf_<nn>_1-95Hz_reg_RASTER.mat      deltaFoF, raster, movements
│   │   └── <date>_GCaMP5_6dpf_<nn>_1-95Hz_reg_ARTIFACTS.mat   not used by the code
│   └── MeCP2/
│       └── <date>_GCaMP5_<n>dpf_MeCP2_null_<nn>_1-95Hz{.mat, _reg_ALL_CELLS.mat, ...}
└── MECP2_VR_Behavior/                              Figure 4 (40 mecp2, 38 wild-type larvae)
    ├── mecp2/
    │   ├── bouts.csv                               one row per bout, all larvae
    │   └── qc.csv                                  larvae excluded by quality control
    └── nacre/
        ├── bouts.csv
        └── qc.csv
```

To use your own data, keep the same folder and file name patterns:

- Files are searched recursively, so subfolders are allowed. The genotype is taken from the
  path, which must contain `WT` or `MeCP2`.
- In Figures 2 and 3, cell data and stimulation files are paired by sorted file name, so each
  larva's files must sort in the same order.
- `bouts.csv` columns used: `file`, `stim`, `category`, `proba`, `sign`, `trial_time`,
  `prey_arc_start_deg`, `foreground_color`. `qc.csv` needs a `file` column.
- The setup render for Supplementary Figure 2A is in the repository:
  `Figure4/3d_render/setup_vr_pc.png`.

## Repository layout

```
mecp2_paper/
├── run_all.py              runs every figure (MATLAB and Python)
├── download_data.py        downloads and verifies the data
├── environment.yml         Python environment (conda-forge)
├── Common/                 shared MATLAB code: bootstrap, decoder, statistics, plotting helpers
├── Figure1/ … Figure3/     F<n>.m (analysis), plot_F<n>.m (layout), panel_*_F<n>.m (panels),
│                           extra/ (exploratory scripts, not used for the paper)
└── Figure4/                F4.py (analysis), plot_F4.py / plot_SupF4.py (layout),
                            panel_* / SupPanel_* (panels), mixed_anova.py, style_F4.py
```