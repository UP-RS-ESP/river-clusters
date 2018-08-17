# create some fake slope area plots to illustrate the clustering

import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
import pandas as pd
import math


def slope_from_drainage_area(areas, ks, theta):
    noise_level = ks/10000
    slopes = [ks*math.pow(x, -theta) for x in areas]
    mean = np.mean(slopes)
    # add some random noise
    noise = np.random.normal(mean, noise_level, size = len(areas))
    #print noise
    slopes=slopes+noise
    #ks * math.pow(areas, theta)
    return slopes

nr_of_profiles = 5
fig, ax = plt.subplots(nrows=nr_of_profiles,ncols=1,sharex=True,sharey=True, figsize =(5,8))
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
colors = ['b', 'r', 'g', 'k', 'orange']

areas = np.linspace(1000,100000,1000)
theta = np.random.uniform(low=0.2, high=0.8, size=5)
ks_values = np.random.uniform(low=0.01, high=1, size=5)

for i,ks in enumerate(ks_values):
    rand_areas = np.random.choice(areas,size=500)
    slopes = slope_from_drainage_area(rand_areas, ks, theta[i])
    ax[i].scatter(rand_areas,slopes, c=colors[i],s=2)
    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    ax[i].set_ylim(0.0001,1)
    #plt.scatter(rand_areas,slopes)

# set axis labels
plt.xlabel('Drainage area (m$^2$)', fontsize=16)
plt.ylabel('Gradient', fontsize=16)
plt.ylim(0.001,0.1)
plt.xscale('log')
plt.yscale('log')
plt.subplots_adjust(left=0.15)
#plt.show()
plt.savefig('/home/clubb/Bitbucket/hrt_workshop/figures/synthetic_sa_plots.png', dpi=300)
plt.clf()

def two_profiles_example():
    """
    Make a SA plot with two different samples as a diagram to show an
    example of the clustering
    """
    
