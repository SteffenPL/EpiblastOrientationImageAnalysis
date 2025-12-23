"""Rotate a composite image into the embryo coordinate frame using rotation metadata."""

import glob
import os
import pathlib
import sys

import imageio.v3 as iio
import numpy as np
import scipy.ndimage as ndi
import skimage as ski
import toml

calibration = np.array([0.1632, 0.1632, 0.1632])
annotations = "ExE_cells"

try:
    snakemake
except NameError:
    snakemake = None

# get our inputs
if snakemake:
    input_image_fn = str(snakemake.input.input_image)
    rotation_metadata_fn = str(snakemake.input.rotation_metadata)
    output_image_fn = str(snakemake.output.rotated_image)
    params = snakemake.params

else:
    embryo_folder = pathlib.Path("/Users/SteffenPlunder/Documents/Workspace/large_datasets/takafumi/epiblast/wildtype/epi90/230914_E525_4-7_nolumen")

    input_image_fn = str(list(embryo_folder.glob("Binned_composite_*.tif"))[0])
    rotation_metadata_fn = str(embryo_folder / "snakemake" / "rotates" / "rotation_metadata.toml")
    output_image_fn = str(embryo_folder / "snakemake" / "rotated" / "rotated_composite_ch1.tif")
    params = {"channel": 1, "order": 1}

    print(f"Process: \n\t{[input_image_fn, rotation_metadata_fn]} \n\t-> {output_image_fn}")

os.makedirs(os.path.dirname(output_image_fn), exist_ok=True)

channel_param = params.get("channel", None)
order_param = int(params.get("order", 1))

metadata = toml.load(rotation_metadata_fn)
R_embryo_to_image = np.asarray(metadata["R_embryo_to_image"], dtype=float)
offset = np.asarray(metadata["offset"], dtype=float)
pixel_size = np.asarray(metadata["pixel_size"], dtype=float)
image_size = np.asarray(metadata["image_size"], dtype=int)

assert np.all(pixel_size[0] == pixel_size), "Non-isotropic pixel size not supported"
factor = int(round(pixel_size[0] / calibration[0]))

R_image_to_embryo = R_embryo_to_image.T  # orthonormal basis => inverse is transpose

# load image
image = iio.imread(input_image_fn)

if image.ndim == 4:
    if channel_param is None:
        raise ValueError("channel parameter is required for 4D images")
    channel_index = int(channel_param)
    if channel_index < 1:
        raise ValueError("channel parameter should be 1-based (>=1)")
    channel_index -= 1  # convert to zero-based
    print(f"Rotating channel {channel_index + 1}")
    image_channel = image[:, channel_index, :, :]
elif image.ndim == 3:
    # Already single-channel 3D; ignore channel parameter
    print("Rotating single-channel 3D image")
    image_channel = image
else:
    raise ValueError(f"Unsupported image ndim={image.ndim}; expected 3 or 4")

# output_shape is stored downscaled in metadata; reconstruct full-resolution shape for transform
new_shape = (image_size * factor).astype(int)

rotated = ndi.affine_transform(
    image_channel,
    R_embryo_to_image,
    offset=offset,
    order=order_param,
    output_shape=new_shape
)

if order_param == 0:
    rotated = (
        ski.transform.rescale(
            rotated,
            1.0 / factor,
            order=0,
            anti_aliasing=False,
            preserve_range=True,
            channel_axis=None,
        ).astype(rotated.dtype)
    )
else:
    rotated = ski.transform.downscale_local_mean(rotated, factor)

iio.imwrite(output_image_fn, rotated)
