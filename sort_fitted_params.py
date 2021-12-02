import pandas as pd
import numpy as np
import sys
import re


def sort_multicol_groups(arr, sortcols, n_groups):
    """
    sort rows of an array in column-groups across, based on the order that
    sorts the numbers in the sortcols

    NOTE: n_groups indicates the TOTAL number of col groups, including the
          sortcols' group!

    NOTE: THE SECOND LINE OF CODE MAKES THE HARD ASSUMPTION THAT THE SORTCOLS
          ARE THE LEFTMOST COL GROUP IN THE ARRAY!
    """
    # get list of the indices needed to sort the sortcols' vals in each row
    sortcols_idxs = np.stack([np.argsort(arr[i,
                                    sortcols]) for i in range(arr.shape[0])])
    # get the list of the indices needed to sort all cols' vals in each row,
    # based on the sortcols
    sort_idxs = np.hstack([sortcols_idxs] + [
                        sortcols_idxs+((i+1)*len(sortcols)) for i in range(n_groups-1)])

    # sort the array using the sort_idxs
    sorted = np.array([arr[np.repeat(i, len(sortcols)*n_groups),
                           sort_idxs[i, :]] for i in range(arr.shape[0])])

    return sorted


infile, outfile = sys.argv[1:3]

df = pd.read_csv(infile, na_filter=False)
cols = [col for col in df.columns if re.search('(^mu)|(^kappa)|(^alpha)', col)]
vals_to_sort = df.loc[:, cols].values
sorted = sort_multicol_groups(vals_to_sort, [0,1,2,3], 3)
df.loc[:, cols] = sorted
df.to_csv(outfile, index=False)
print('%s written' % outfile)
