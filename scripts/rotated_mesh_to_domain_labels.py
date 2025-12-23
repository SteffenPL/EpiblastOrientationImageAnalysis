"""Create a label image from the rotated embryo mesh assigning domain IDs.

Labels:
    0 -> background
    1 -> inside embryo mesh below cut-plane height
    2 -> inside embryo mesh above cut-plane height
"""

from pathlib import Path
import numpy as np
import toml
import trimesh
import imageio.v3 as iio


try:
    snakemake
except NameError:
    snakemake = None

if snakemake:
    mesh_fn = str(snakemake.input.mesh)
    metadata_fn = str(snakemake.input.metadata)
    output_label_fn = str(snakemake.output.label_image)
else:
    sample_dir = Path("/Users/SteffenPlunder/Documents/Workspace/large_datasets/takafumi/epiblast/wildtype/epi65/230601_E500_1-2")
    mesh_fn = str(sample_dir / "snakemake" / "rotated" / "rotated_embryo_mesh.stl")
    metadata_fn = str(sample_dir / "snakemake" / "rotated" / "rotation_metadata.toml")
    output_label_fn = str(sample_dir / "snakemake" / "rotated" / "rotated_domain_labels.tif")

metadata = toml.load(metadata_fn)
image_size = np.asarray(metadata["image_size"], dtype=int)
cut_plane_height = float(metadata["rotated_cut_plane_height"])

mesh = trimesh.load(mesh_fn)

# Voxelize the mesh at pitch=1 (mesh is already in rotated voxel coordinates)
voxel_grid = mesh.voxelized(pitch=1.0).fill()

# convert voxel centers to integer indices in the rotated image grid
centers = voxel_grid.points  # (N, 3) world coords of occupied voxel centers
indices = np.round(centers).astype(int)

# keep only indices inside the target volume
in_bounds = (
    (indices[:, 0] >= 0)
    & (indices[:, 1] >= 0)
    & (indices[:, 2] >= 0)
    & (indices[:, 0] < image_size[0])
    & (indices[:, 1] < image_size[1])
    & (indices[:, 2] < image_size[2])
)
indices = indices[in_bounds]

labels = np.zeros(image_size, dtype=np.uint8)

# build cut-plane masks on the full image grid
y_coords = np.arange(image_size[1], dtype=float)[None, :, None]
below_mask = y_coords < cut_plane_height
above_mask = ~below_mask

if indices.size > 0:
    inside_mask = np.zeros(image_size, dtype=bool)
    inside_mask[indices[:, 0], indices[:, 1], indices[:, 2]] = True
    labels[inside_mask & below_mask] = 1
    labels[inside_mask & above_mask] = 2

iio.imwrite(output_label_fn, labels)
