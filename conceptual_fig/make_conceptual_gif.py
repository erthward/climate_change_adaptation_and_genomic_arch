import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.stats import norm
from itertools import cycle

plt.rc('animation', html='html5')

save = False

class AnimatedPlot(object):
    """An animated plot using matplotlib.animations.FuncAnimation.
       Stolen and tweaked from:
           https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot
    """
    def __init__(self, time_series, stable, density, n_individs):
        self.time_series_list = time_series
        self.time_series = cycle(time_series)
        self.stable = stable
        self.density = density
        self.curr_timestep = next(self.time_series)
        self.n_individs = n_individs

        self.stream = self.data_stream()

        # Setup the figure and axes...
        self.fig = plt.figure()
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        self.ax1.set_title('shifting gradient')
        self.ax2.set_title('stable gradient')
        self.ax1.axes.xaxis.set_visible(False)
        self.ax1.axes.yaxis.set_visible(False)
        self.ax2.axes.xaxis.set_visible(False)
        self.ax2.axes.yaxis.set_visible(False)

        # holders for the artists
        self.im1 = None
        self.im2 = None
        self.scat1 = None
        self.scat2 = None

        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=500, 
                                          init_func=self.setup_plot, blit=True)

    def setup_plot(self):
        """Initial drawing of the plot."""
        x, y = next(self.stream)
        c = x + np.random.normal(0, 1, len(x))
        self.im1 = self.ax1.imshow(self.curr_timestep, vmin=0, vmax=1, cmap='coolwarm')
        self.im2 = self.ax2.imshow(self.stable, vmin=0, vmax=1, cmap='PiYG')
        self.scat1 = self.ax1.scatter(x, y, c=c, cmap="coolwarm_r", edgecolor="k")
        self.scat2 = self.ax2.scatter(x, y, c=x, cmap='PiYG_r', edgecolor='k')
        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.im1, self.im2, self.scat1, self.scat2

    def data_stream(self):
        """Draw individuals' locations."""
        x, y = draw_individs(self.density[0,:], self.n_individs)
        while True:
            x, y = draw_individs(self.density[0,:], self.n_individs)
            yield x, y

    def update(self, i):
        """Update the plot."""


        x, y = next(self.stream)

        # Set x and y data...
        self.scat1.set_offsets(np.vstack((x, y)).T)
        self.scat2.set_offsets(np.vstack((x, y)).T)

        # Update the raster
        self.curr_timestep = next(self.time_series)
        self.im1.set_data(self.curr_timestep)

        # Set colors..
        new_ax1_colors = np.array([1-self.curr_timestep[0,:][int(val)] for val in x]) * 50
        self.scat1.set_array(new_ax1_colors)
        self.scat2.set_array(x)

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.im1, self.scat1, self.scat2


# create a time series of changing rasters (corresponding to first trait)
time_series = [np.vstack([np.hstack((np.ones(i),
                                     np.linspace(1, 0, 25),
                                     np.zeros(25-i))) for _ in range(25)])
               for i in range(25)]

# create a stable raster (corresponding to second trait)
stable = np.vstack([np.linspace(1, 0, 50) for _ in range(25)])

# create a stable density raster (centered, normal distribution)
dens_vals = [norm.pdf(n) for n in np.linspace(-2.5, 2.5, 50)]
dens_vals = [(val - min(dens_vals))/(
              max(dens_vals) - min(dens_vals)) for val in dens_vals]
density = np.vstack([dens_vals for _ in range(25)])

# create a stack of shifting density rasters, as an alternative scenario
dens_stack = [(ts_elem + stable)/2 for ts_elem in time_series]


# function to draw random individual values
def draw_individs(density_vals, n_individs):
    individ_x_counts = np.array([np.random.binomial(n_individs, p,
                                                1)[0] for p in density_vals])
    individ_xs = np.hstack([np.random.normal(i, 0.5, 
                            n) for i, n in enumerate(individ_x_counts)]) - 0.5
    individ_ys = np.random.uniform(low=0, high=25, size=len(individ_xs)) - 0.5
    individ_xs = np.clip(individ_xs, 0.0001, 49.999)
    individ_ys = np.clip(individ_ys, 0.0001, 49.999)
    return (individ_xs, individ_ys)

fig = AnimatedPlot(time_series, stable, density, 20)
plt.show()
if save:
    fig.ani.save('conceptual.gif', writer='imagemagick', fps=10)

"""

# make a gif of the stable-density scenario

fig_stable_dens = plt.figure('stable density')
ax1 = fig_stable_dens.add_subplot(211)
ax1.set_title('shifting gradient')
ax2 = fig_stable_dens.add_subplot(212)
ax2.set_title('stable gradient')
ax1.imshow(time_series[0], cmap='coolwarm', vmin=0, vmax=1)
ax1.xticks=([])
ax1.yticks=([])
ax2.imshow(stable, cmap='PiYG', vmin=0, vmax=1)
ax1.axes.xaxis.set_visible(False)
ax1.axes.yaxis.set_visible(False)
ax2.axes.xaxis.set_visible(False)
ax2.axes.yaxis.set_visible(False)

def init():
    ax1.imshow(time_series[0], cmap='coolwarm', vmin=0, vmax=1)
    ax1.xticks=([])
    ax1.yticks=([])
    return fig_stable_dens,

ax_dict = {}

def animate(i):
    print('Now making frame %i...' % (i+1))
    if i == 0:
        xs, ys = draw_individs(dens_vals)
        ax_dict[0] = ax1.scatter(xs, ys, c=xs + np.random.normal(0, 1, len(xs)),
                             cmap='coolwarm_r', edgecolors='black', s=0)
        ax_dict[1] = ax2.scatter(xs, ys, c=xs + np.random.normal(0, 1, len(xs)),
                             cmap='PiYG_r', edgecolors='black', s=0)

    ax1.imshow(time_series[i], cmap='coolwarm', vmin=0, vmax=1)
    ax1.xticks=([])
    ax1.yticks=([])
    xs, ys = draw_individs(dens_vals)
    new_ax1_colors = np.array([1-time_series[i][0,:][int(x)] for x in xs]) * 50
    ax_dict[0].set_offsets(np.vstack((xs, ys)).T)
    ax_dict[1].set_offsets(np.vstack((xs, ys)).T)
    ax_dict[0].set_array(new_ax1_colors)
    ax_dict[1].set_array(xs)
    ax_dict[0].set_sizes(np.array([25 for _ in range(len(xs))]))
    ax_dict[1].set_sizes(np.array([25 for _ in range(len(xs))]))
    return fig_stable_dens,

anim = animation.FuncAnimation(fig_stable_dens, animate, init_func=init,
                               frames=25, interval=25*20, blit=True)
anim.save('TEST.gif', writer='imagemagick', fps=60)


# make a gif of the shifting-density scenario
"""
