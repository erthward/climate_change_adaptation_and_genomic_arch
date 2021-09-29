import numpy as np
import matplotlib.pyplot as plt

scape_size = 100

ys = xs = np.arange(0, scape_size, 1)
dirs = []
origs = []
dests = []
#for rand_origin in range(200):
#    # draw a random origin
#    rand_orig = np.random.uniform(0, scape_size, 2)
#    origs.append(rand_orig)
#    for x in xs:
#        for y in ys:
#            # get the effective coords of the movement relative to the origin
#            eff_x = x - rand_orig[0]
#            eff_y = y - rand_orig[1]
#            dir = np.arctan2(eff_y, eff_x)/(2*np.pi)*360
#            if dir<0:
#                dir+=360
#            dirs.append(dir)
for rand_origin in range(20000):
    # draw a random origin
    rand_orig = np.random.uniform(0, scape_size, 2)
    rand_dest = np.random.uniform(0, scape_size, 2)
    # allow them to move from and to points with the population's bulk density
    # on the landscape instead?...
    #rand_orig = np.random.uniform(20, scape_size-20, 2)
    #rand_dest = np.random.uniform(40, scape_size, 2)
    origs.append(rand_orig)
    dests.append(rand_dest)
    # get the effective coords of the movement relative to the origin
    eff_x = rand_dest[0] - rand_orig[0]
    eff_y = rand_dest[1] - rand_orig[1]
    dir = np.arctan2(eff_y, eff_x)/(2*np.pi)*360
    if dir<0:
        dir+=360
    dirs.append(dir)


plt.hist(dirs, bins = 30)
plt.show()
