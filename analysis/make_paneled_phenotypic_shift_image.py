import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import os, re

# loosen PIL Image-size limit, to prevent DecompressionBombWarning error
#Image.MAX_IMAGE_PIXELS = None

# load all images
#analysis_dir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev'
analysis_dir = '/global/scratch/users/drewhart/ch2/output/analysis'
im_files = [f for f in os.listdir(analysis_dir) if re.search(
                                                    '^phenotypic_shift.*', f)]
ims = {f:Image.open(f, 'r') for f in im_files}

# get size of an image
width, height = [*ims.values()][0].size

# set sizes of margins to be cropped and replaced with white
crop_left = 250
crop_top = 175
crop_right = -350
crop_bott = -185


# generate the paneled output image
nrow=3
ncol=3
out_im = Image.new('RGB',(int(width*ncol*1.01), int(height*nrow*1.01)))

# function to pad an image
# adapted from: https://note.nkmk.me/en/python-pillow-add-margin-expand-canvas/
def add_margin(pil_img, left, top, right, bottom, color='white'):
    width, height = pil_img.size
    new_width = width + right + left
    new_height = height + top + bottom
    result = Image.new(pil_img.mode, (new_width, new_height), color)
    result.paste(pil_img, (left, top))
    return result


# fill the output image with the input images
for i, linkage in enumerate(['independent', 'weak', 'strong']):
    for j, genicity in enumerate([4, 20, 100]):

        # get the current image
        curr_im_filename_patt = '^phenotypic_shift_L%s_G\d*%i\.png' % (linkage,
                                                                       genicity)
        curr_im = [v for k, v in ims.items() if re.search(curr_im_filename_patt,
                                                          k)]
        assert len(curr_im) == 1
        curr_im = curr_im[0]

        # crop the current image's margins and replace them with white,
        # to avoid repeating mariginal info in output image
        crop_box = [0,0,curr_im.size[0], curr_im.size[1]]
        margins_to_add = [0,0,0,0]
        # crop the lefts, if not in the left row
        if j > 0:
            crop_box[0] = crop_left
            margins_to_add[0] = crop_left
        # crop the image tops, if not in the top row
        #curr_im = curr_im.crop((0, 100, 0, 0))
        if i > 0:
            crop_box[1] = crop_top
            margins_to_add[1] = crop_left
        # crop the rights, if not in the right row
        if j < 2:
            crop_box[2] = crop_box[2] + crop_right
            margins_to_add[2] = -crop_right
        # crop the bottoms, if not in the bottom rows
        if i < 2:
            crop_box[3] = crop_box[3] + crop_bott
            margins_to_add[3] = -crop_bott

        curr_im = add_margin(curr_im.crop(crop_box), *margins_to_add)

        # add some space between

        px, py = j*width, i*height
        out_im.paste(curr_im,(px,py))

# save image to file
out_im.save(os.path.join(analysis_dir, 'phenotypic_shift_grid_fig.jpg'))
