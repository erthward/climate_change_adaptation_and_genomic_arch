import tskit
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import bisect


# for a sample of nodes and a sample of loci, get the sequentially ordered
# dict of lineage dicts (containing each node in a child node's lineage as the
# key and its birth-time and birth-location as a length-2 list of values
def get_lineage_dicts(tc, nodes, loci, drop_before_sim=True):
    # make sure the TableCollection is sorted, then get the TreeSequence
    try:
        ts = tc.tree_sequence()
    except Exception:
        raise Exception(("The species' TableCollection must be sorted and "
                         "simplified before this method is called."))
    # get the tree numbers corresponding to each locus
    treenums = get_treenums(ts, loci)
    # get the trees
    trees = ts.aslist()
    # create the output data structure
    lineage_dicts = {}
    # get the lineage_dict for each locus
    for loc, treenum in zip(loci, treenums):
        # store the lineage dicts for the current locus
        lineage_dicts_curr_loc = get_lineage_dicts_one_tree(tc, trees[treenum],
                                    nodes, drop_before_sim=drop_before_sim)
        lineage_dicts[loc] = lineage_dicts_curr_loc
    return lineage_dicts



# for a sample of nodes, and for a given tree,
# get dicts of each node's lineage nodes with their birth locs and birth times
def get_lineage_dicts_one_tree(tc, tree, nodes, drop_before_sim=True):
    #create an output master dict
    lineage_dicts = {}
    # get each sample node's lineage dict
    for node in nodes:
        lineage = get_lineage(tree, [node])
        lineage_dict = dict(zip(lineage, get_lineage_times_and_locs(tc,
                                                                    lineage)))
        # optionally drop all nodes from before the simulation (i.e. which
        # were added to TableCollection by backward-time simulation with ms)
        if drop_before_sim:
            lineage_dict = {node: val for node, val in lineage_dict.items(
                                                            ) if val[0] < 0}
        lineage_dicts[node] = lineage_dict
    return lineage_dicts


# recursive algorithm for getting all parents of a node
def get_lineage(tree, nodes):
    parent = tree.parent(nodes[-1])
    if parent == tskit.NULL:
        return(nodes)
    else:
        nodes.append(parent)
        return(get_lineage(tree, nodes))


# get birth-locations and times of nodes along a pedigree
def get_lineage_times_and_locs(tc, lineage):
    locs = []
    times = []
    for node in lineage:
        time = tc.nodes[node].time
        times.append(time)
        individ = tc.nodes[node].individual
        loc = tc.individuals[individ].location
        locs.append(loc)
    return(zip(times, locs))


# get a dict of the interval (i.e. tree) number
# keyed to each locus in a list of loci
def get_treenums(ts, loci, as_dict=False):
    interval_right_ends = [t.get_interval()[1] for t in ts.trees()]
    treenums = [bisect.bisect_left(interval_right_ends, loc) for loc in loci]
    if as_dict:
        treenums = {loc: treenum for loc, treenum in zip(loci, treenums)}
    return(treenums)


# plot lineage for a given node and locus
def plot_lineages(spp, locus, land, nodes=None, individs=None, phenotype=None,
                  lyr_num=0, jitter=True):
    assert nodes is not None or individs is not None, ("Either the nodes "
                "argument or the individs argument must be provided. "
                "If both are provided, only the nodes argument's value "
                "will be used.")
    # sort and simplify the TableCollection, then grab it and its TreeSequence
    tc = spp._tc
    try:
        ts = tc.tree_sequence()
    except Exception:
        raise Exception(("The species' TableCollection must be sorted and "
                         "simplified before this method is called."))
    # get the individs' nodes, if necessary
    if nodes is None:
        nodes = np.hstack([[*spp[ind]._nodes_tab_ids.values(
                                                    )] for ind in individs])
    # get the tree for this locus
    tree = ts.aslist()[get_treenums(ts, [locus])[0]]
    # get the lineage_dict (with birth times and birth locations)
    lin_dict = get_lineage_dicts_one_tree(tc, tree, nodes)
    # create start-color values for nodes' separate lineage tracks
    start_cols = [mpl.colors.to_hex(plt.cm.Set1(n)) for n in np.linspace(0, 
                                                                    0.85, 8)]
    #start_cols = ['#ffa114', '#ff1956', '#e419ff', '#4f6fff', '#4fffdf',
    #             '#87ff4f']
    # plot the species, either with or without phenotype-painting
    if phenotype is None:
        spp._plot(lyr_num=lyr_num, land=land)
    else:
        spp._plot_phenotype(phenotype, land=land)
    ax = plt.gca()
    # extract and plot the series of points for each node
    for i, node in enumerate(nodes):
        start_col = start_cols[i % len(start_cols)]
        locs = np.vstack([v[1] for v in lin_dict[node].values()])
        if jitter:
            locs = locs + np.random.normal(0, 0.01,
                                           size=locs.size).reshape(locs.shape)
        # create linear interpolation of colors for plotting
        color_nums = np.int8(np.linspace(0, 100, locs.shape[0]-1))
        colors = [color_variant(start_col, num) for num in color_nums]
        # create a linear interpolation of linewidths
        linewidths = np.linspace(3, 0.85, locs.shape[0]-1)
        for n, color in enumerate(colors):
            ax.plot(locs[n:n+2, 0], locs[n:n+2, 1], linestyle='--', marker='o',
                    markersize=25**(1/2), color=color, linewidth=linewidths[n])
    plt.show()


# return an arbitrarily lighter version of a color
# (stolen and tweaked from Chase Seibert: https://chase-seibert.github.io/blog/
# 2011/07/29/python-calculate-lighterdarker-rgb-colors.html)
def color_variant(hex_color, brightness_offset=1):
    """ takes a color like #87c95f and produces a lighter or darker variant
    """
    if len(hex_color) != 7:
        raise Exception(("Passed %s into color_variant(), needs to be "
                        "in #87c95f format.") % hex_color)
    rgb_hex = [hex_color[x:x+2] for x in [1, 3, 5]]
    new_rgb_int = [int(hex_value,
                       16) + brightness_offset for hex_value in rgb_hex]
    # make sure new values are between 0 and 255
    new_rgb_int = [min([255, max([0, i])]) for i in new_rgb_int]
    # hex() produces "0x88", we want just "88"
    new_hex_int = [hex(i)[2:] for i in new_rgb_int]
    new_hex_int = [str(i).zfill(2) for i in new_hex_int]
    new_hex = "#" + "".join(new_hex_int)
    return new_hex


# calculate the angular direction of gene flow along a locus' lineage
def calc_gene_flow_dir_dist_time(lineage):
    # get x and y distances between beginning and ending points
    beg_loc = [*lineage.values()][0][1]
    end_loc = [*lineage.values()][-1][1]
    x_diff, y_diff = [end_loc[i] - beg_loc[i] for i in range(2)]
    # get the counterclockwise angle, expressed in degrees
    # from the vector (X,Y) = (1,0),
    # with 0 to 180 in quadrants 1 & 2, 0 to -180 in quadrants 3 & 4
    ang = np.rad2deg(np.arctan2(y_diff, x_diff))
    # convert to all positive values, 0 - 360
    if ang < 0:
        ang += 360
    #convert to all positive clockwise angles, expressed as 0 at compass north
    # (i.e. angles clockwise from vector (X,Y) = (0,1)
    ang = (-ang + 90) % 360

    dist = np.sqrt(x_diff**2 + y_diff**2)

    beg_time = [*lineage.values()][0][0]
    end_time = [*lineage.values()][-1][0]
    time = end_time - beg_time

    return ang, dist, time

