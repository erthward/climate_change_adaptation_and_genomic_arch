#import rpy2.robjects as robjects
#from rpy2.robjects import pandas2ri
from rpy2.rinterface_lib import openrlib

import scipy
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

with openrlib.rlock:
    import rpy2_subprocess as rp2sp
    rp2sp.main()


