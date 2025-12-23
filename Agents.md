# Agents guide

- Purpose: Snakemake workflow for embryo/epiblast microscopy that builds meshes, rotates data, and outputs labels/metadata.
- Core tools: Snakemake orchestrates Python scripts using numpy/scipy, image I/O, and trimesh for mesh handling. Python >=3.12 (`pyproject.toml` lists full deps).

## Key files
- `workflow/Snakefile`: Picks the first existing path in `inputdirs` from `snakemake_config.yaml`, discovers samples via `Selected_*.tif`, and skips any with `"epi15"` in the folder name. The default `all` target builds meshes, rotation metadata, rotated meshes, rotated labels, rotated cut/tip metadata, and domain labels for each sample (channels 1–4 rotates are available but commented out). Creates a `combined_rotation_metadata.toml` by merging per-sample metadata.
- `scripts/*.py`: Rule actions (meshing, rotations, metadata, labels). See script notes below.

## Script notes (scripts/)
- `label_to_embryo_mesh.py`: Marching cubes on `Selected_*.tif` labels → `embryo_mesh.stl` (assumes isotropic 0.1632 µm voxels).
- `label_to_cell_meshes.py`: Builds per-label glTF meshes (colored) plus embryo mesh and cut_plane point cloud; step_size=4 marching cubes.
- `embryo_mesh_to_coordinates.py`: Derives rotation/origin/cut_plane_height from embryo mesh + `ExE_cells/cut_plane.csv`, `tip_point.csv`, `front_points.csv`; auto tip-point saved to `auto_tip_point.csv` unless `ignore_auto_tip_point.csv` exists.
- `rotate_metadata.py`: Projects cut_plane/tip into rotated frame; writes `rotated_cut_plane.csv`, `rotated_tip_point.csv`, and `rotation_metadata.toml` (pixel_size, image_size, origins, offset, R matrices). Assumes downscale factor=2 and uses `Binned_composite_*.tif` for shape.
- `rotate_image_into_coordinates.py`: Re-samples composite/label images into embryo coordinates using affine; supports 3D or 4D with `channel` param; downscales by 2 (rescale for labels).
- `rotate_mesh_into_coordinates.py`: Applies rotation metadata to embryo mesh and saves rotated STL.
- `rotated_mesh_to_domain_labels.py`: Voxelizes rotated mesh; labels 1 = inside below cut-plane height, 2 = inside above, 0 = background.
- `rotate_pts_into_coordinates copy.py`: Transforms xyz columns of a CSV into embryo frame; Snakemake references `rotate_points_into_coordinates.py` which is missing—likely needs this script renamed.
- `back_rotate_pts_from_coordinates.py`: Inverse transform of rotated points back to image coordinates using rotation metadata.
- `downscale_image.py`: Utility to downscale TIF via `downscale_local_mean` (default factor=2).
- `compress_repo.py`: Zips repo contents (tracked + unignored) to `tmp/snapshot_<repo>_<date>.zip`.

## Running
- Install deps with `uv sync` (or your env manager of choice), then run: `snakemake -s workflow/Snakefile --cores 4 -p`.
- Without an active venv, use: `uv run snakemake -s workflow/Snakefile --cores 4 -p`.
- Dry run: add `-n` (or `-nkp` to also print commands). Override data roots inline if needed: `--config inputdirs='["/path/one","/path/two"]'`.

## Notes
- Outputs live under each sample’s `snakemake/` folder; `.gitignore` excludes these. Calibration hardcoded to 0.1632 µm; rotated outputs assume downscale factor=2.
- Inputs per sample: `Selected_*.tif`, optionally `Binned_composite_*.tif`, and `ExE_cells/cut_plane.csv`, `tip_point.csv`, `front_points.csv` (+ optional `ignore_auto_tip_point.csv`). Samples containing `"epi15"` are excluded by default.
- No automated tests; sanity-check on a small sample dataset.
