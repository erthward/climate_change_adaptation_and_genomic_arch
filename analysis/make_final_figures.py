import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import os, re, sys


# function to pad an image
# adapted from: https://note.nkmk.me/en/python-pillow-add-margin-expand-canvas/
def add_margin(pil_img, left, top, right, bottom, color='white'):
    width, height = pil_img.size
    new_width = width + right + left
    new_height = height + top + bottom
    result = Image.new(pil_img.mode, (new_width, new_height), color)
    result.paste(pil_img, (left, top))
    return result


# loosen PIL Image-size limit, to prevent DecompressionBombWarning error
Image.MAX_IMAGE_PIXELS = None

# make paired Nt and mean fit time series and boxplot figures
for stat in ['mean_fit', 'Nt']:
    ts_files = ['ch2_%s_loREDUND_over_time.jpg' % stat,
                'ch2_%s_hiREDUND_over_time.jpg' % stat,
               ]
    box_files = ['boxplot_delta_%s_loREDUND.jpg' % re.sub('mean_', '', stat),
                 'boxplot_delta_%s_hiREDUND.jpg' % re.sub('mean_', '', stat),
                ]
    ts_ims = {f:Image.open(f, 'r') for f in ts_files}
    box_ims = {f:Image.open(f, 'r') for f in box_files}

    # get sizes of images
    ts_width, ts_height = [*ts_ims.values()][0].size
    box_width, box_height = [*box_ims.values()][0].size

    # set sizes of margins to be cropped and replaced with white
    crop_left = 0
    crop_top = 0
    crop_right = 0
    crop_bott = 0


    # generate the paneled output image
    nrow=2
    ncol=2
    out_im = Image.new('RGB', ts_width*1.1 + box_width*1.1)


    # fill the output image with the input images
    for i, redundancy in ['lo', 'hi']:

                out_im.paste(curr_im,(px,py))

# save image to file
if plot_type == 'DENS':
    out_im.save(os.path.join(analysis_dir,
                    'pop_density_shift_grid_fig_%sREDUND.jpg' % redundancy))
else:
    out_im.save(os.path.join(analysis_dir,
        'phenotypic_shift_grid_fig%s_%sREDUND.jpg' % (
                                            '_NULL' * (nullness == 'null'),
                                            redundancy)))
