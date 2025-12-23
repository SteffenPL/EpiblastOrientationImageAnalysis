"""Downscale a TIF image by an integer factor (for quick previews or lighter processing)."""

from pathlib import Path

import skimage as ski
from imageio.v3 import imread, imwrite

try:
    snakemake
except NameError:
    snakemake = None

if snakemake:
    input_image_fn = str(snakemake.input.image)
    output_image_fn = str(snakemake.output.downscaled)
    factor = snakemake.params.get("factor", 2)

else:
    fn = Path("/Users/SteffenPlunder/Documents/Workspace/large_datasets/takafumi/epiblast/wildtype/epi90/230914_E525_4-7_nolumen/snakemake/rotated/rotated_composite_ch1.tif")
    input_image_fn = str(fn)
    output_image_fn = fn.parent / ("binned_" + fn.name)
    output_image_fn = str(output_image_fn)

    factor = 2

img = imread(input_image_fn)
img_downscaled = ski.transform.downscale_local_mean(img, factor)
imwrite(output_image_fn, img_downscaled)
