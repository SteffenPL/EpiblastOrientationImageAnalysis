# EpiblastOrientationImageAnalysis

Repository analyzing epiblast orientation maps supporting

- "Boundary-guided cell alignment drives mouse epiblast maturation." by T. Ichikawa, P. C. Guruciaga, S. Hu, S. Plunder, M. Makino, M. Hamaji, A. Stokkermans, S. Yoshida, T. Hiiragi, A. Erzberger.
- "Boundary geometry controls a topological defect transition that determines lumen nucleation in embryonic development" by P. C. Guruciaga, T. Ichikawa, S. Plunder, T. Hiiragi, A. Erzberger.

## Description

Snakemake workflow for processing embryo/epiblast microscopy data and generating meshes, rotated metadata, and domain labels.

## Setup with uv
- Install `uv` if needed: `curl -LsSf https://astral.sh/uv/install.sh | sh`.
- Create/sync the environment from `pyproject.toml`: `uv sync`.
- Activate the env when working (optional if you use `uv run`): `source .venv/bin/activate` (or `.\.venv\Scripts\activate` on Windows).

## Running Snakemake
- Provide one or more data paths in `snakemake_config.yaml` under `inputdirs` (first existing path is used). You can also override on the command line.
- Run the workflow: `snakemake -s workflow/Snakefile --cores 4 -p`.
- If the venv is not active, use uv to run it: `uv run snakemake -s workflow/Snakefile --cores 4 -p`.
- Dry run (no file changes): append `-n` for a plan, or `-nkp` to also show commands without touching files.
- Add/override data paths inline when needed: append `--config inputdirs='["/path/one", "/path/two"]'`.
- Outputs appear inside each sample folder under `snakemake/` (meshes, rotated metadata, domain labels, rotated composites when enabled).

## Required input files
- The workflow reads `snakemake_config.yaml` and picks the first `inputdirs` entry that exists on disk; every discovered sample lives two directories under that root (e.g., `<inputdir>/<batch>/<sample>/`). Any directory whose path contains `epi15` is ignored automatically.
- Each sample directory **must contain exactly one** `Selected_*.tif` label stack (3D). This volume is used for both embryo and per-cell meshing—multiple matches will cause Snakemake to abort.
- Each sample also needs **at least one** `Binned_composite_*.tif` image. The first alphabetical match is used to resample up to four composite channels via `rotate_image_into_coordinates_ch[1-4]` as `(z, channel, y, x)`.
- Provide manual annotations inside `ExE_cells/` within every sample:
  - `cut_plane.csv`: ≥3 points describing the cut plane as exported from Imaris/napari (columns `index,axis-0,axis-1,axis-2` in voxel coordinates).
  - `tip_point.csv`: one or more rows marking the embryo tip; multiple rows will be averaged, then optionally replaced by an automatically detected tip.
  - `front_points.csv`: anterior reference points; averaged to orient the embryo.
  - Optional `ignore_auto_tip_point.csv`: when present, `embryo_mesh_to_coordinates.py` keeps the manual tip and writes the auto-detected tip into this file; otherwise the auto tip replaces the manual one and is saved to `auto_tip_point.csv`.
- All CSV coordinates are interpreted in pixel units and internally converted to microns using the hard-coded 0.1632 µm isotropic calibration, so ensure the annotation files match the voxel spacing of the `Selected_*` and `Binned_composite_*` images.

# Post-processing

To generate orientation map plots, we refer to the notebooks in `/submodules/EpiblastOrientationMaps` or (https://github.com/SteffenPL/EpiblastOrientationMaps) and the readme therein.