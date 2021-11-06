import scipy
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

# example of fitting distributional params to a hist of data,
# then running KS test
# (NOTE: technically invalidates KS test assumptions, bc params should not
#  have been estimated from the data; consider Anderson-Darling or Cramer-von
#  Mises tests instead)
# from: https://medium.com/@amirarsalan.rajabi/distribution-fitting-with-
#       python-scipy-bb70a42c0aed
sample = np.random.normal(5,2,1000)
dist = stats.norm
params = dist.fit(sample)
print(params)
# p >> 0.05, so H0 that data are ~N(5,2) is not rejected
stats.kstest(sample, "norm", params)



# read data from some of my simulation output
data = pd.read_csv('output_PID-182509_DIR_short.csv',
                   na_filter=False)

non = [n for n in data[data.genicity == 100][data.linkage=='independent'][data.neutrality=='nonneut'][data.nullness=='non-null'][data.it==0]['dir']]
nons = []
for x in non:
    try:
        nons.append(float(x)/360*(2*np.pi))
    except Exception as e:
        print(x)
        print(e)

# fit von Mises mixture distribution
# (adapting code from:
    #https://stackoverflow.com/questions/35990467/fit-mixture-of-two-gaussian
    #-normal-distributions-to-a-histogram-from-one-set-of#_=_
y,x,_=hist(nons,100,alpha=.3,label='non-null')
x=(x[1:]+x[:-1])/2 # for len(x)==len(y)

def vonmises(x,mu,kappa):
    # NOTE: scipy.special.k0 is the modified Bessel fn of order 0
    return (np.e**(kappa*np.cos(x-mu)))/(2*np.pi*scipy.special.k0(kappa))

def vonmises_4mix(x,mu1,kappa1,mu2,kappa2,mu3,kappa3,mu4,kappa4):
    out = (vonmises(x,mu1,kappa1)+
           vonmises(x,mu2,kappa2)+
           vonmises(x,mu3,kappa3)+
           vonmises(x,mu4,kappa4))
    return out


expected=(1,0,1,np.pi/2,1,np.pi,1,3*np.pi/2)
params,cov=curve_fit(vonmises_4mix,x,y,expected)
sigma=sqrt(diag(cov))
plot(x,vonmises_4mix(x,*params),color='red',lw=3,label='model')
legend()
print(params,'\n',sigma)


