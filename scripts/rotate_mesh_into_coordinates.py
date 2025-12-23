"""Rotate the embryo mesh into the coordinate frame defined by the rotation metadata."""

from pathlib import Path

import numpy as np
import toml
import trimesh

try:
    snakemake
except NameError:
    snakemake = None

if snakemake:
    embryo_fn = str(snakemake.input.embryo_mesh)
    rotation_data_fn = str(snakemake.input.rotation_data)
    output_label_fn = str(snakemake.output.label_image)

else:
    sample_dir = Path("/Users/SteffenPlunder/Documents/Workspace/large_datasets/takafumi/epiblast/wildtype/epi65/230601_E500_1-2")

    # inputs
    embryo_fn = str(sample_dir / "snakemake" / "embryo_mesh.stl")
    rotation_data_fn = str(sample_dir / "snakemake" / "rotated" / "rotation_metadata.toml")

    # outputs
    rotated_embryo_fn = str(sample_dir / "snakemake" / "rotated" / "rotated_embryo_mesh.stl")
    output_label_fn = rotated_embryo_fn


# load embryo mesh
embryo_mesh = trimesh.load(embryo_fn)


# load rotation data
rot_data = toml.load(rotation_data_fn)
R_image_to_embryo = np.asarray(rot_data["R_image_to_embryo"], dtype=float)
offset = np.asarray(rot_data["offset"], dtype=float)
pixel_size = np.asarray(rot_data["pixel_size"], dtype=float)

assert np.all(pixel_size[0] == pixel_size), "Non-isotropic pixel size not supported"
factor = pixel_size[0] / 0.1632  #from microns to original voxel size

# apply the same rotation/offset/downscale as used for the rotated images

rotated_vertices = (R_image_to_embryo @ (embryo_mesh.vertices.T/0.1632 - offset[:, np.newaxis])).T
rotated_vertices /= factor
embryo_mesh.vertices = rotated_vertices.astype(np.float32)

# save scaled embryo mesh
embryo_mesh.export(output_label_fn)
