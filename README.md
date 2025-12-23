# ti-epiblast-cell-organization

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
