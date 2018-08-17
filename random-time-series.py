# Plot some random time series for a presentation

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

# number of time series
n = 5

# random time series of random length
ts = []

fig, ax = plt.subplots(nrows=n,ncols=1,sharex=True,sharey=True, figsize =(5,8))
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

x = np.linspace(0,100,1)

colors = ['b', 'r', 'g', 'k', 'orange']

for i in range(n):
    # random length rlen
    y = np.random.random(100)
    ax[i].plot(y, c = colors[i])

plt.xlabel('Time (e.g. hours)', fontsize=16)
plt.ylabel('Some kind of metric', fontsize=16)
plt.savefig('/home/clubb/Bitbucket/hrt_workshop/figures/random_ts.png', dpi=300, transparent=True)
