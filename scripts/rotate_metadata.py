"""Rotate cut plane and tip point metadata into the embryo coordinate frame using rotation_metadata.toml."""

import os
import pathlib

import numpy as np
import pandas as pd
import toml

calibration = np.array([0.1632, 0.1632, 0.1632])

try:
    snakemake
except NameError:
    snakemake = None

if snakemake:
    rotation_metadata_fn = str(snakemake.input.rotation_metadata)
    cut_plane_fn = str(snakemake.input.cut_plane)
    tip_point_fn = str(snakemake.input.tip_point)
    rotated_cut_plane_fn = str(snakemake.output.rotated_cut_plane)
    rotated_tip_point_fn = str(snakemake.output.rotated_tip_point)
    params = snakemake.params
else:
    embryo_folder = pathlib.Path("/Users/SteffenPlunder/Documents/Workspace/large_datasets/takafumi/epiblast/wildtype/epi90/230914_E525_4-7_nolumen")
    rotation_metadata_fn = str(embryo_folder / "snakemake" / "rotated" / "rotation_metadata.toml")
    cut_plane_fn = str(embryo_folder / "ExE_cells" / "cut_plane.csv")
    tip_point_fn = str(embryo_folder / "ExE_cells" / "tip_point.csv")
    rotated_cut_plane_fn = str(embryo_folder / "snakemake" / "rotated" / "rotated_cut_plane.csv")
    rotated_tip_point_fn = str(embryo_folder / "snakemake" / "rotated" / "rotated_tip_point.csv")
    params = {}
    print(
        f"Process metadata: \n\t{[rotation_metadata_fn, cut_plane_fn, tip_point_fn]} "
        f"\n\t-> {[rotated_cut_plane_fn, rotated_tip_point_fn]}"
    )

os.makedirs(os.path.dirname(rotated_cut_plane_fn), exist_ok=True)

metadata = toml.load(rotation_metadata_fn)
R_image_to_embryo = np.asarray(metadata["R_image_to_embryo"], dtype=float)
offset = np.asarray(metadata["offset"], dtype=float)
pixel_size = np.asarray(metadata["pixel_size"], dtype=float)

assert np.all(pixel_size[0] == pixel_size), "Non-isotropic pixel size not supported"
factor = pixel_size[0] / calibration[0]

cut_plane = pd.read_csv(cut_plane_fn)
cut_plane_points = cut_plane[["axis-0", "axis-1", "axis-2"]].to_numpy()
rotated_cut_plane_points = (R_image_to_embryo @ (cut_plane_points.T - offset[:, np.newaxis])).T

rotated_cut_plane_df = pd.DataFrame(
    data={
        "index": np.arange(rotated_cut_plane_points.shape[0]),
        "axis-0": rotated_cut_plane_points[:, 0] / factor,
        "axis-1": rotated_cut_plane_points[:, 1] / factor,
        "axis-2": rotated_cut_plane_points[:, 2] / factor,
    }
)
rotated_cut_plane_df.to_csv(rotated_cut_plane_fn, index=False)

tip_points = pd.read_csv(tip_point_fn)
tip_points_array = tip_points[["axis-0", "axis-1", "axis-2"]].to_numpy()
rotated_tip_points = (R_image_to_embryo @ (tip_points_array.T - offset[:, np.newaxis])).T

rotated_tip_points_df = pd.DataFrame(
    data={
        "index": np.arange(rotated_tip_points.shape[0]),
        "axis-0": rotated_tip_points[:, 0] / factor,
        "axis-1": rotated_tip_points[:, 1] / factor,
        "axis-2": rotated_tip_points[:, 2] / factor,
    }
)
rotated_tip_points_df.to_csv(rotated_tip_point_fn, index=False)
