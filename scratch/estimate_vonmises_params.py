import numpy as np

def calc_mean_std_angs(angs):
    """
    Calculates the mean and std of a list or 1d array of angles (in radians).

    Does this by decomposing each angle into its component
    unit vectors, calculating the means and root mean squared deviations
    of the x and y component vectors, then composing the two mean vectors
    into a mean angle and the two root mean squared vectors into an angular
    standard deviation.
    """
    # get component vectors
    xs = [np.cos(ang) for ang in angs]
    ys = [np.sin(ang) for ang in angs]
    # get mean component vectors
    mean_x = np.mean(xs)
    mean_y = np.mean(ys)
    # get root mean squared deviation component vectors
    rmsdev_x = np.sqrt(np.mean([(x - mean_x)**2 for x in xs]))
    rmsdev_y = np.sqrt(np.mean([(y - mean_y)**2 for y in ys]))
    print(rmsdev_x, rmsdev_y)
    # compose mean vectors
    mean_ang = np.arctan(mean_y/mean_x)
    std_ang = np.arctan(rmsdev_y/rmsdev_x)

    print(mean_y, mean_x)
    std_ang = np.sqrt(-np.log(mean_y**2 + mean_x*2))
    return (mean_ang, std_ang)


def estimate_vonmises_params(angs, p=2):
    angs = [*np.array(angs)[~np.isnan(angs)]]
    xi = np.stack((np.cos(angs), np.sin(angs)))
    xbar = np.sum(xi, axis=1)/len(angs)
    Rbar = np.sqrt(np.sum([val**2 for val in xbar]))
    muhat = xbar/Rbar
    kappahat = (Rbar*(p-Rbar**2))/(1-Rbar**2)
    d = 1 - np.sum([muhat.dot(xi[:,i]) for i in range(xi.shape[1])])/len(angs) 
    sigmahat = np.sqrt((d)/(len(angs)*(Rbar**2)))
    muhat_ang = np.arctan2(*muhat[::-1])
    return muhat_ang, kappahat, sigmahat
