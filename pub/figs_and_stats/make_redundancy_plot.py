import numpy as np
import matplotlib.pyplot as plt
from math import comb

plt.close('all')

save_it = True


def calc_redundancy(genicity, log=True):
    """
    Laruson et al. Box 1, eqxn ii, implemented for diploid life cycle with 0.5
    (rather than 0) serving as the baseline phenotype.
    """
    if genicity in [4, 20, 100]:
        redundancy = 'lo'
    elif genicity in [8, 40, 200]:
        redundancy = 'hi'
    else:
        raise ValueError
    α = 1/genicity
    if redundancy == 'hi':
        α *= 2
    k = (2*genicity)+1
    base_genotype = 0.5
    min_z = base_genotype - (genicity/2*α)
    max_z = base_genotype + (genicity/2*α)
    zs = np.arange(min_z, max_z+(α/2), α/2)
    assert len(zs) == k
    assert np.allclose(np.mean(zs), 0.5)
    n_genos = [comb(2*genicity, val) for val in range(k)]
    if log:
        n_genos = np.log(np.array(n_genos).astype(float))
    return zs, n_genos


def plot_redundancy(ax, genicity, color, log=True):
    ax.plot(*calc_redundancy(genicity, log=log), color=color)
    for pheno_thresh in range(2):
        ax.axvline(x=pheno_thresh, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1],
                   color='red', linewidth=0.5, linestyle='-.', alpha=0.5)


colors = {'lo': '#32e3ba', # lighter aqua
          'hi': '#0a6e57', # darker aqua
         }
fig = plt.figure(figsize=(5,11))
for i, genicity in enumerate([4, 20, 100]):
    ax = fig.add_subplot(3, 1, i+1)
    ax.set_title('%i loci (low redundancy) or\n%i loci (high redundancy)' %
                 (genicity, genicity*2), fontdict={'fontsize':17})
    for redundancy in ['lo', 'hi']:
        # double genicity if high-redundancy
        if redundancy == 'hi':
            genicity = genicity * 2
        # plot log of number of genotypes per phenotype
        plot_redundancy(ax, genicity, color=colors[redundancy])
    if i == 2:
        ax.set_xlabel('phenotype', fontdict={'fontsize':16})
    else:
        ax.set_xlabel('')
    ax.set_ylabel('log(# of genotypes)', fontdict={'fontsize':16})
    ax.tick_params(labelsize=12)
    # label number of loci
    ax.text(0.83, 0.33*ax.get_ylim()[1], genicity, color=colors['lo'], size=10)
    ax.text(1.2, 0.6*ax.get_ylim()[1], genicity*2, color=colors['hi'], size=10)


fig.subplots_adjust(right=0.94, left=0.15, bottom=0.08, top=0.94, hspace=0.65)
fig.show()

if save_it:
    fig.savefig('redundancy_fig.png', dpi=700)

