"""Derive rotation metadata directly from the embryo mesh, cut plane, tip point, and composite image."""

import glob
import os
import pathlib
import sys

import imageio.v3 as iio
import numpy as np
import pandas as pd
import toml
import trimesh

calibration = np.array([0.1632, 0.1632, 0.1632])
annotations = "ExE_cells"

try:
    snakemake
except NameError:
    snakemake = None

if snakemake:
    mesh_fn = str(snakemake.input.embryo_mesh)
    cut_plane_fn = str(snakemake.input.cut_plane)
    tip_point_fn = str(snakemake.input.tip_point)
    binned_composite_fn = str(snakemake.input.binned_composite)
    rotation_metadata_fn = str(snakemake.output.rotation_metadata)
    params = snakemake.params
else:
    embryo_folder = pathlib.Path("/Users/SteffenPlunder/Documents/Workspace/large_datasets/takafumi/epiblast/wildtype/epi90/230913_E525_3-3")

    mesh_fn = str(embryo_folder / "snakemake" / "embryo_mesh.stl")
    cut_plane_fn = str(embryo_folder / "ExE_cells" / "cut_plane.csv")
    tip_point_fn = str(embryo_folder / "ExE_cells" / "tip_point.csv")
    binned_composite_fn = str(list(embryo_folder.glob("Binned_composite_*.tif"))[0])
    rotation_metadata_fn = str(embryo_folder / "snakemake" / "rotated" / "rotation_metadata.toml")
    params = {}

    print(
        f"Process: \n\t{[mesh_fn, cut_plane_fn, tip_point_fn, binned_composite_fn]} "
        f"\n\t-> {rotation_metadata_fn}"
    )

inputdir = str(pathlib.Path(mesh_fn).parent.parent)


def unit_vector(vector, fallback=0.0):
    """Return unit vector, or fallback-filled vector if length is zero."""
    l = np.linalg.norm(vector)
    if l > 0:
        vector = vector / l
    else:
        vector = np.array([fallback] * 3)
    return vector


def angle_between(v1, v2):
    """Returns the angle in radians between vectors 'v1' and 'v2' (in deg)."""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return 180 / np.pi * np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def toHomMat(mat):
    return np.block([[mat, np.zeros((mat.shape[0], 1))], [np.zeros(mat.shape[1]), 1]])


def compute_auto_tip_point(epi_mesh, cut_plane_pts, tip_pt):
    # Get cutting plane
    cut_origin, cut_normal = trimesh.points.plane_fit(cut_plane_pts)

    # Flip normal vector of plane if necessary
    if np.dot(tip_pt - cut_origin, cut_normal) > 0:
        cut_normal = -cut_normal

    # Find point with maximal distance
    dists = np.dot(epi_mesh.vertices, -cut_normal)

    return epi_mesh.vertices[np.argmax(dists)]


def compute_embryo_coordinates(epi_mesh, cut_plane_pts, front_pts, tip_pt, extra_info=None):
    # Step 2: Get cutting plane
    cut_origin, cut_normal = trimesh.points.plane_fit(cut_plane_pts)

    # Step 2.b: Get embryo bottom
    bottom_pt = tip_pt
    front_pts = np.mean(front_pts, axis=0)

    # Step 2.c: Flip normal vector of plane if necessary
    if np.dot(bottom_pt - cut_origin, cut_normal) > 0:
        cut_normal = -cut_normal

    # Step 3: Move cutton plane a bit down
    upper_plane_origin = cut_origin - cut_normal * 10

    # Step 4: Intersect mesh with cutting plane and get the intersection line
    upper_slice = epi_mesh.convex_hull.section(plane_origin=upper_plane_origin, plane_normal=cut_normal)

    # Step 5: Get the centroid of the upper slice
    upper_centroid = upper_slice.centroid
    upper_centroid = upper_centroid + cut_normal * 10  # move it back up to the cut plane

    tip_points = trimesh.proximity.closest_point(epi_mesh, bottom_pt.reshape(1, 3))
    mesh_tip = tip_points[0][0]

    # Step 6: Create coordinates
    origin = epi_mesh.centroid
    up_arrow = origin - mesh_tip  # up_arrow = cut_normal
    back_arrow = upper_centroid - front_pts.mean(axis=0)

    # make orthogonal
    right_arrow = -np.cross(back_arrow, up_arrow)
    back_arrow = -np.cross(up_arrow, right_arrow)

    # normalize
    up_arrow = up_arrow / np.linalg.norm(up_arrow)
    back_arrow = back_arrow / np.linalg.norm(back_arrow)
    right_arrow = right_arrow / np.linalg.norm(right_arrow)
    down_arrow = -up_arrow

    # Step 7: Create a coordinate system
    rotation = toHomMat(np.array([back_arrow, down_arrow, right_arrow]).T)

    # Step 8: We use the mesh centroid as the coordinate system origin, and keep track where the cutting plane is.

    cut_plane_height = np.dot(upper_centroid - origin, down_arrow)

    if extra_info is not None:
        correction_angle = angle_between(cut_normal, up_arrow)
        extra_info["correction_angle"] = correction_angle

    return rotation, origin, cut_plane_height


epi_mesh = trimesh.load(mesh_fn)

cut_plane_pts = np.loadtxt(os.path.join(inputdir, annotations, "cut_plane.csv"), skiprows=1, delimiter=",")[:, 1:]
cut_plane_pts *= calibration

tp = np.loadtxt(os.path.join(inputdir, annotations, "tip_point.csv"), skiprows=1, delimiter=",")

if tp.ndim == 1:
    tip_pt = tp[1:]
else:
    tip_pt = np.mean(tp[:, 1:], axis=0)

tip_pt *= calibration

auto_tp = compute_auto_tip_point(epi_mesh, cut_plane_pts, tip_pt)

# write tip point to file
tp_dataframe = pd.DataFrame(
    data={
        "index": [0.0],
        "axis-0": [auto_tp[0] / calibration[0]],
        "axis-1": [auto_tp[1] / calibration[1]],
        "axis-2": [auto_tp[2] / calibration[2]],
    }
)

# decided if we should ignore the automatic tip point for this given mesh or not:
ignore_automic_tp = os.path.exists(os.path.join(inputdir, annotations, "ignore_auto_tip_point.csv"))
if ignore_automic_tp:
    # store the new tip point in the ignored file
    auto_tp_fn = os.path.join(inputdir, annotations, "ignore_auto_tip_point.csv")
    tp_dataframe.to_csv(auto_tp_fn)
else:
    # store the new tip point in the auto tip point file
    auto_tp_fn = os.path.join(inputdir, annotations, "auto_tip_point.csv")
    tp_dataframe.to_csv(auto_tp_fn)

    # use automatically generate tip point
    tip_point_difference = np.linalg.norm(tip_pt - auto_tp)
    tip_pt = auto_tp

front_pts = np.loadtxt(os.path.join(inputdir, annotations, "front_points.csv"), skiprows=1, delimiter=",")[:, 1:]
front_pts *= calibration

rotation, origin, cut_plane_height = compute_embryo_coordinates(epi_mesh, cut_plane_pts, front_pts, tip_pt)

assert np.unique(calibration).size == 1  # only uniform calibration supported right now

rotation = rotation[:3, :3]
origin = origin / calibration[0]
cut_plane_height = cut_plane_height / calibration[0]

# load composite to derive rotated volume shape and offset
image = iio.imread(binned_composite_fn)
if image.ndim == 4:
    image_channel = image[:, 0, :, :]
elif image.ndim == 3:
    image_channel = image
else:
    raise ValueError(f"Unsupported image ndim={image.ndim}; expected 3 or 4")

R_embryo_to_image = rotation
R_image_to_embryo = R_embryo_to_image.T

xm, ym, zm = image_channel.shape
corners = np.array(
    [
        [0, 0, 0],
        [xm, 0, 0],
        [0, ym, 0],
        [0, 0, zm],
        [xm, ym, 0],
        [xm, 0, zm],
        [0, ym, zm],
        [xm, ym, zm],
    ],
    dtype=float,
)

corners_embryo = np.array((R_image_to_embryo @ (corners).T))
min_embryo = corners_embryo.min(axis=1)
max_embryo = corners_embryo.max(axis=1)
new_shape = np.ceil(max_embryo - min_embryo).astype(int).flatten()
offset = np.array(R_embryo_to_image @ min_embryo).flatten()

factor = 2  # current downscale factor; pixel size scales by this
scaled_calibration = calibration * factor
rotated_shape = np.ceil(new_shape / factor).astype(int)
new_origin = R_image_to_embryo @ (origin.T - offset.T)
scaled_origin = new_origin / factor
new_cut_plane_height = (new_origin[1] + cut_plane_height) / factor

metadata = {
    "pixel_size": scaled_calibration.tolist(),
    "image_size": rotated_shape.tolist(),
    "origin_image": origin.tolist(),
    "origin_rotated": scaled_origin.tolist(),
    "offset": offset.tolist(),
    "rotated_cut_plane_height": float(new_cut_plane_height),
    "R_image_to_embryo": R_image_to_embryo.tolist(),
    "R_embryo_to_image": R_embryo_to_image.tolist(),
}

os.makedirs(os.path.dirname(rotation_metadata_fn), exist_ok=True)
with open(rotation_metadata_fn, "w") as f:
    toml.dump(metadata, f)
