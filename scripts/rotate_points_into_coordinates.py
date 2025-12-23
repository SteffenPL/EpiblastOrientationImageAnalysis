"""Rotate the embryo mesh into the coordinate frame defined by the rotation metadata."""

from pathlib import Path

import numpy as np
import pandas as pd
import toml

try:
    snakemake
except NameError:
    snakemake = None

if snakemake:
    points_fn = str(snakemake.input.points)
    rotation_data_fn = str(snakemake.input.rotation_data)
    rotated_points_fn = str(snakemake.output.rotated_points)

else:
    sample_dir = Path("/Users/SteffenPlunder/Documents/Workspace/large_datasets/takafumi/epiblast/wildtype/epi65/230601_E500_1-2")

    # in
    points_fn = str(sample_dir / "ExE_cells" / "cut_plane.csv")
    rotation_data_fn = str(sample_dir / "snakemake" / "rotated" / "rotation_metadata.toml")
    # out 
    rotated_points_fn = str(sample_dir / "snakemake" / "rotated" / "cut_plane.csv")


# load pts (last three columns are x,y,z)
pts_df = pd.read_csv(points_fn)
pts_in = pts_df.iloc[:, -3:].to_numpy()

# load rotation data
rot_data = toml.load(rotation_data_fn)
R_image_to_embryo = np.asarray(rot_data["R_image_to_embryo"], dtype=float)
offset = np.asarray(rot_data["offset"], dtype=float)
pixel_size = np.asarray(rot_data["pixel_size"], dtype=float)
calibration = 0.1632 

assert np.all(pixel_size[0] == pixel_size), "Non-isotropic pixel size not supported"

def rotate_points(pts, R, offset, calibration, pixel_size):
    """Rotate points into embryo coordinate frame."""
    factor = pixel_size[0] / calibration  # scaling between original voxel size and downscaled/rotated image voxel size

    rotated_pts = (R @ (pts.T - offset[:, np.newaxis])).T
    rotated_pts /= factor
    return rotated_pts

rotated_pts = rotate_points(pts_in, R_image_to_embryo, offset, calibration, pixel_size)

# makedirs if needed
output_dir = Path(rotated_points_fn).parent
output_dir.mkdir(parents=True, exist_ok=True)

print(f"Saving rotated points to {rotated_points_fn}")

# save rotated points while preserving other columns
pts_df.iloc[:, -3:] = rotated_pts
pts_df.to_csv(rotated_points_fn, index=False)
