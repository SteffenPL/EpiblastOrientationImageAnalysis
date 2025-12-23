"""Convert labeled embryo voxels into a single embryo surface mesh."""

import glob
import os
import pathlib
import sys

import numpy as np
import tifffile
import trimesh

calibration = np.array([0.1632, 0.1632, 0.1632])

try:
    snakemake
except NameError:
    snakemake = None

if snakemake:
    label_image_fn = str(snakemake.input.label_image)
    output_mesh_fn = str(snakemake.output.embryo_mesh)
    params = snakemake.params

    assert output_mesh_fn.endswith(".stl")

else:
    try:
        dir = pathlib.Path("/Users/SteffenPlunder/Documents/Workspace/large_datasets/takafumi/epiblast/231016_ltgb1_embryos/")
        embryo_folder = list(dir.glob("*"))[0]
        label_image_fn = list(embryo_folder.glob("Selected_*.tif"))[0]
        output_mesh_fn = embryo_folder / "snakemake" / "embryo_mesh.stl"
        params = {"voxel_size": 1.0}

    except:
        label_image_fn = ""
        output_mesh_fn = ""
        params = {}

labelimg = tifffile.imread(label_image_fn)
mesh = trimesh.voxel.ops.matrix_to_marching_cubes(labelimg, pitch=calibration)

# move center of mass to origin
# mesh.vertices -= mesh.center_mass

os.makedirs(os.path.dirname(output_mesh_fn), exist_ok=True)
mesh.export(output_mesh_fn)
