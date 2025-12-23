"""Convert labeled voxel data into individual cell meshes and save as glTF."""

import glob
import os
import pathlib
import sys

import distinctipy
import numpy as np
import pandas as pd
import tifffile
import trimesh
from skimage.measure import marching_cubes
from trimesh.exchange import gltf

calibration = np.array([0.1632, 0.1632, 0.1632])

try:
    snakemake
except NameError:
    snakemake = None

if snakemake:
    label_image_fn = str(snakemake.input.label_image)
    output_mesh_fn = str(snakemake.output.cell_meshes)
    params = snakemake.params

    assert output_mesh_fn.endswith(".gltf")

else:
    
    dir = pathlib.Path("/Users/SteffenPlunder/Documents/Workspace/large_datasets/takafumi/epiblast/wildtype/epi90/")
    embryo_folder = list(dir.glob("*_E*"))[0]

    label_image_fn = list(embryo_folder.glob("Selected_*.tif"))[0]
    output_mesh_fn = embryo_folder / "snakemake" / "cell_meshes" / "cell_meshes.gltf"
    params = {"voxel_size": 1.0}

    print(f"Process: \n\t{label_image_fn} \n\t-> {output_mesh_fn}")

labelimg = tifffile.imread(label_image_fn)
labels = np.unique(labelimg)
labels = labels[labels != 0]  # remove background and border

labels = labels
imgsize = labelimg.shape
offset = np.array(imgsize) / 2.0

scene = trimesh.Scene()
colors = distinctipy.get_colors(len(labels))


for k in range(len(labels)):
    l = labels[k]
    color = colors[k]
    print(f"Processing: {k} / {len(labels)}")

    labelimg_single = labelimg == l

    verts, faces, _, _ = marching_cubes(labelimg_single, level=0.5, allow_degenerate=False, step_size=4, spacing=calibration)
    mesh = trimesh.Trimesh(vertices=verts, faces=faces)

    mesh.vertices -= offset

    # export 
    mesh.visual.vertex_colors = [int(c*255) for c in color] + [255]
    scene.add_geometry(mesh, geom_name="label" + str(l))

verts, faces, _, _ = marching_cubes(labelimg > 0, level=0.5, allow_degenerate=False, step_size=4, spacing=calibration)
embryo_mesh = trimesh.Trimesh(vertices=verts, faces=faces)
embryo_mesh.vertices -= offset

scene.add_geometry(embryo_mesh, geom_name="embryo")

embryo_dir = os.path.dirname(label_image_fn)
cut_plane = os.path.join(embryo_dir, "ExE_cells", "cut_plane.csv")
cut_plane_df = pd.read_csv(cut_plane)

# save as points in the scene
points = cut_plane_df[["axis-0", "axis-1", "axis-2"]].to_numpy()
points -= offset
point_cloud = trimesh.points.PointCloud(points, colors=[255, 0, 0, 255])
scene.add_geometry(point_cloud, geom_name="cut_plane")


# save mesh
os.makedirs(os.path.dirname(output_mesh_fn), exist_ok=True)
scene.export(output_mesh_fn)
