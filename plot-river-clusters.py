#---------------------------------------------------------------------#
# Clustering of river profiles
# Developed by Fiona Clubb
#              Bodo Bookhagen
#              Aljoscha Rheinwalt
# University of Potsdam
#---------------------------------------------------------------------#

# setting backend to run on server
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from glob import glob
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
from scipy import stats
from scipy.ndimage.interpolation import shift
from CorrCoef import Pearson
import math
import LSDPlottingTools as LSDP
from LSDMapFigure.PlottingRaster import MapFigure
import sys
from collections import defaultdict
import os
from matplotlib import ticker

# Set up fonts for plots
label_size = 10
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = label_size

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to do some river profile clustering for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("You also need the -fname flag which will give the prefix of the raster files.")
    print("For help type:")
    print("   python plot_river_clusters.py -h\n")
    print("=======================================================================\n\n ")

#---------------------------------------------------------------------#
# ANALYSIS FUNCTIONS
#---------------------------------------------------------------------#
def find_nearest_idx(array,value):
    """
    Given a value, find the index of the point in the array which is closest
    to that value.
    Author: FJC
    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def find_difference_between_arrays(x, y):
    """
    Function to calculate a difference between two np arrays with the same
    dimensions. If values are all positive, scales from 0 to 1: 1 = no difference, 0 = bad.
    1 - (||(X-Y)/(X+Y)|| / sqrt(n))
    """
    n = len(x)
    #print n
    #print x, y
    num = x - y
    #print num
    den = x + y
    div = np.divide(num,den)
    norm = np.linalg.norm(div)
    #mean_diff = abs(np.mean(x - y))
    #print norm
    #print np.sqrt(n)
    diff = 1-(norm/np.sqrt(n))
    #diff = 1 - mean_diff
    #print diff
    return diff

def ProfilesRegularDistance(df, profile_len = 1000, step=2, slope_window_size=25):
    """
    This function takes the dataframe of river profiles and creates an array
    that can be used for the time series clustering.  For each profile, the slope
    is assigned to a common distance step, in metres (default=1m)

    Args:
        df: pandas dataframe from the river profile csv.
        profile_len (int): number of data points in each profile (= distance from source)
        step (int): step size that you want in metres. This should be greater than the maximum
        possible spacing between your points (max_spacing = sqrt(2 * DataRes^2)). Default = 2
        slope_window_size (int): window over which slope was calculated

    Returns: array of size (n_profiles, profile_len) containing the slopes along each profile

    Author: FJC
    """
    print("The profile length is {} m, any smaller channels will be removed.".format(profile_len))
    print("Assigning the profiles a common distance step of {} m".format(step))
    # find the slope radius. We don't calculate slope for the first or last few nodes, which
    # are below the min required window size. Therefore reg dist needs to be smaller
    slope_radius = ((slope_window_size-1)/2)

    # create new array of regularly spaced differences
    reg_dist = np.arange(step, profile_len, step)


    # find the minimum length that the array can be (profile length/root2)
    min_length = profile_len/(math.sqrt(2))

    # loop through the dataframe and store the data for each profile as an array of
    # slopes and distances
    profiles = []
    source_ids = df['id'].unique()
    final_sources = []
    for i, source in enumerate(source_ids):
        this_df = df[df['id'] == source]
        this_df = this_df[np.isnan(this_df['slope']) == False]  # remove nans
        if not this_df.empty:
            slopes = this_df['slope'].as_matrix()[::-1] #need to reverse these as the distances need to be sorted
            distances = this_df['distance_from_outlet'].as_matrix()[::-1]
            if (distances.max() >= profile_len):
                profiles.append((distances, slopes))
                #thinned_df = thinned_df.append(this_df)
                final_sources.append(source)

    # now create the 2d array to store the data
    n_profiles = len(profiles)
    data = np.empty((n_profiles, len(reg_dist)))

    # create a new dataframe for storing the data about the selected profiles
    thinned_df = pd.DataFrame()

    # make a plot of the gradient vs. distance from source aligned to the outlet with
    # the resampled distance frame
    # set up a figure
    # fig = plt.figure(1, facecolor='white')
    # gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    # ax = fig.add_subplot(gs[5:100,10:95])
    #
    # loop through the profiles. For each point in the regularly spaced array,
    # find the index of the closest point in the distance array. Then use this to
    # assign the slope to the regularly spaced array
    for i, p in enumerate(profiles):
        reg_slope = []
        for d in reg_dist:
            idx = find_nearest_idx(p[0], d)
            reg_slope.append(p[1][idx])
            this_dist = p[0][idx]
            # get this distance and append the regular distance to the thinned df
            thinned_df = thinned_df.append(df.loc[(df['id']==final_sources[i]) & (df['distance_from_outlet'] == p[0][idx])])
        thinned_df.loc[(thinned_df['id'] == final_sources[i]), 'reg_dist'] = reg_dist
        data[i] = reg_slope

        # plot this profile
        # ax.plot(reg_dist[::-1], reg_slope, lw=1)

    # write the thinned_df to output in case we want to reload
    thinned_df.to_csv(DataDirectory+fname_prefix+'_profiles_upstream_reg_dist.csv',index=False)
    #
    # # now save the figure
    # ax.set_xlabel('Distance from source (m)')
    # ax.set_ylabel('Gradient')
    # plt.savefig(DataDirectory+fname_prefix+'_profiles_upstream_reg_dist.png', dpi=300)
    # plt.clf()

    return thinned_df, data

def ProfilesRegDistVaryingLength(df, profile_len=4,step=2, slope_window_size=25):
    """
    This function takes the dataframe of river profiles and creates an array
    that can be used for the time series clustering.  For each profile, the slope
    is assigned to a common distance step, in metres (default=1m).
    This can be used on profiles of varying distance: instead of creating the array here
    we just write the regularly spaced distances to the dataframe.

    Args:
        df: pandas dataframe from the river profile csv.
        profile_len (int): number of unique points you need to keep a profile
        step (int): step size that you want in metres. This should be greater than the maximum
        possible spacing between your points (max_spacing = sqrt(2 * DataRes^2)). Default = 2
        slope_window_size (int): window over which slope was calculated

    Returns: df with regularly spaced distances.

    Author: FJC
    """
    print("Assigning the profiles a common distance step of {} m".format(step))
    # find the slope radius. We don't calculate slope for the first or last few nodes, which
    # are below the min required window size. Therefore reg dist needs to be smaller
    slope_radius = ((slope_window_size-1)/2)

    # make a plot of the gradient vs. distance from source aligned to the outlet with
    # the resampled distance frame
    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # loop through the dataframe and store the data for each profile as an array of
    # slopes and distances
    source_ids = df['id'].unique()
    print source_ids
    rows_list = []
    for i, source in enumerate(source_ids):
        this_df = df[df['id'] == source]
        # create new array of regularly spaced differences
        reg_dist = np.arange(step, int(this_df.distance_from_outlet.max())+step, step)

        if not this_df.empty:
            df_array = this_df.as_matrix()[::-1]
            slopes = df_array[:,-1]
            distances = df_array[:,5]
            reg_slope = []
            for d in reg_dist:
                # find the index of the nearest point in the distance array
                idx = find_nearest_idx(distances, d)
                reg_slope.append(slopes[idx])
                # find the row that corresponds to this distance
                this_row = df_array[idx,:]
                this_row = np.append(this_row,d)

                rows_list.append(this_row)

            # plot this profile
            ax.plot(reg_dist, reg_slope, lw=1)

    # change the array back to a dataframe
    cols = list(df.columns.values)
    cols.append('reg_dist')
    thinned_df = pd.DataFrame(data=rows_list,columns=cols)

    # now remove non-unique profiles
    thinned_df = RemoveNonUniqueProfiles(thinned_df)
    thinned_df = RemoveProfilesWithShortUniqueSection(thinned_df, profile_len)

    # write the thinned_df to output in case we want to reload
    thinned_df.to_csv(DataDirectory+fname_prefix+'_profiles_upstream_reg_dist_var_length.csv', index=False)

    # now save the figure
    ax.set_xlabel('Distance from outlet (m)')
    ax.set_ylabel('Gradient')
    plt.savefig(DataDirectory+fname_prefix+'_profiles_upstream_reg_dist_var_length.png', dpi=300)
    plt.clf()

    return thinned_df

def ClusterProfiles(df, profile_len=100, step=2, min_corr=0.5, method='complete'):
    """
    Cluster the profiles based on gradient and distance from source.
    Aggolmerative clustering, see here for more info:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage

    Args:
        df: pandas dataframe from the river profile csv.
        profile_len (int): number of data points in each profile (= distance from source)
        step (int): the spacing in metres between the data points that you want, default = 1m
        min_corr (float): minimum correlation threshold for clustering
        method (str): clustering method to use, see scipy docs. Can be 'single', 'complete', 'average',
        'weighted', 'centroid', 'median', or 'ward'. Default is 'complete'.

    Author: AR, FJC
    """
    print ("Now I'm going to do some hierarchical clustering...")

    # get the data from the dataframe into the right format for clustering
    sources = df['id'].unique()
    #print sources
    pts_in_profile = len(df[df['id'] == sources[0]])
    data = np.empty((len(sources), pts_in_profile))

    for i, src in enumerate(sources):
        this_df = df[df['id'] == src]
        data[i] = this_df['slope']
    #print data


    # we could have a look at the ranks too ..
    # correlations
    cc = Pearson(data)

    # distances
    dd = np.arccos(cc)

    # do agglomerative clustering by stepwise pair matching
    # based on angle between scalar products of time series
    ln = linkage(dd, method=method)

    # define threshold for cluster determination
    thr = np.arccos(min_corr)

    # compute cluster indices
    cl = fcluster(ln, thr, criterion = 'distance')
    print("I've finished! I found {} clusters for you :)".format(cl.max()))
    print cl

    set_link_color_palette(colors)

    source_ids = df['id'].unique()

    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    R = dendrogram(ln, color_threshold=thr, above_threshold_color=threshold_color)

    plt.axhline(y = thr, color = 'r', ls = '--')
    plt.savefig(DataDirectory+fname_prefix+"_upstream_dendrogram.png", dpi=300)
    plt.clf()

    for i,id in enumerate(source_ids):
        df.loc[df.id==id, 'cluster_id'] = cl[i]

    return df

def ClusterProfilesVaryingLength(df, profile_len=100, step=2, min_corr=0.5, method='ward'):
    """
    Cluster the profiles based on gradient and distance from source. This works for profiles of varying length.
    Aggolmerative clustering, see here for more info:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage

    Args:
        df: pandas dataframe from the river profile csv.
        min_corr (float): minimum correlation threshold for clustering
        method (str): clustering method to use, see scipy docs. Can be 'single', 'complete', 'average',
        'weighted', 'centroid', 'median', or 'ward'. Default is 'ward'.

    Author: AR, FJC
    """
    print ("Now I'm going to do some hierarchical clustering...")
    np.set_printoptions(threshold='nan')

    #sort the dataframe based on max distance from outlet for each source id.

    # get the data from the dataframe into the right format for clustering
    sources = df['id'].unique()
    n = len(sources)
    print sources
    print n
    data = []

    for i, src in enumerate(sources):
        this_df = df[df['id'] == src]
        data.append(this_df['slope'].as_matrix())

    # correlation coefficients
    cc = np.zeros(int(n * (n - 1) / 2))
    #cc = []
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            tsi = data[i]
            tsj = data[j]
            if len(tsi) > len(tsj):
                tsi = tsi[len(tsi)-len(tsj):]  # Fiona - testing just clustering the top of each profile
                #tsi = tsi[:len(tsj)]
            else:
                #tsj = tsj[:len(tsi)]
                tsj = tsj[len(tsj)-len(tsi):]
            dts = tsi - tsj
            l = 0
            while dts[l] == 0:
                l += 1
            # # try to just cluster over the top 25% of each of the profiles
            #new_l = int(0.25*l)+l
            tsi, tsj = tsi[l:], tsj[l:]
            #tsi, tsj = tsi[new_l:], tsj[new_l:]
            cc[k] = np.corrcoef(tsi, tsj)[0, 1]
            k += 1

    #print np.isnan(cc)
    #print cc
    # distances
    dd = np.arccos(cc)
    #print dd
    print len(dd)
    # do agglomerative clustering by stepwise pair matching
    # based on angle between scalar products of time series
    ln = linkage(dd, method=method)

    # make a plot of the distance vs number of clusters. Use this to determine
    # the threshold
    thr = PlotDistanceVsNClusters(ln)

    # define threshold for cluster determination
    #thr = np.arccos(min_corr)
    #thr = 1.7

    # compute cluster indices
    cl = fcluster(ln, thr, criterion = 'distance')
    print("I've finished! I found {} clusters for you :)".format(cl.max()))
    print cl

    set_link_color_palette(colors)

    source_ids = df['id'].unique()

    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    R = dendrogram(ln, color_threshold=thr+0.00001, above_threshold_color=threshold_color)

    plt.axhline(y = thr, color = 'r', ls = '--')
    plt.savefig(DataDirectory+fname_prefix+"_upstream_dendrogram.png", dpi=300)
    plt.clf()

    for i,id in enumerate(source_ids):
        df.loc[df.id==id, 'cluster_id'] = cl[i]

    return df

def ClusterProfilesDrainageArea(df, profile_len=100, step=2, method='ward'):
    """
    Cluster the profiles based on gradient and drainage area.
    Aggolmerative clustering, see here for more info:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage

    Args:
        df: pandas dataframe from the river profile csv.
        method (str): clustering method to use, see scipy docs. Can be 'single', 'complete', 'average',
        'weighted', 'centroid', 'median', or 'ward'. Default is 'ward'.

    Author: FJC, AR
    """
    print ("Now I'm going to do some hierarchical clustering by drainage area...")
    #np.set_printoptions(threshold='nan')

    # get the data from the dataframe into the right format for clustering
    sources = df['id'].unique()
    n = len(sources)

    all_areas = df['drainage_area'].as_matrix()
    #all_areas = np.log10(all_areas)
    #sort the areas
    sorted_areas = np.sort(all_areas)
    print (len(sorted_areas))
    reg_areas = np.unique(sorted_areas[0:-1:30])
    print (len(reg_areas))

    # create matrix for the data
    data = np.empty((len(sources), len(reg_areas)))

    # now for each profile, find the nearest point on the trunk areas array and assign the slope value to this.
    for x, src in enumerate(sources):
        this_df = df[df['id'] == src]
        reg_slopes = np.empty(len(reg_areas))
        reg_slopes[:] = np.nan

        #df_array = this_df.as_matrix()[::-1]
        slopes = this_df['slope'].as_matrix()
        areas = this_df['drainage_area'].as_matrix()
        #slopes = np.log10(slopes)
        #areas = np.log10(areas)
        for idx, i in enumerate(areas):
            #print i
            idx_j = find_nearest_idx(reg_areas, i)
            #print ('THIS AREA: ', trunk_areas[idx_j])
            reg_slopes[idx_j] = slopes[idx]
        data[x] = reg_slopes

    # correlation coefficients
    cc = np.zeros(int(n * (n - 1) / 2))

    k = 0
    for i in range(n):
        for j in range(i+1, n):
            tsi = data[i]
            tsj = data[j]
            new_tsi = []
            new_tsj = []
            # print ("LEN OLD SERIES 1:", np.count_nonzero(~np.isnan(tsi)))
            # print ("LEN OLD SERIES 2:", np.count_nonzero(~np.isnan(tsj)))
            # remove any areas where there isn't data in both time series
            for x in range(len(tsi)):
                if not (np.isnan(tsi[x])) and not (np.isnan(tsj[x])):
                    new_tsi.append(tsi[x])
                    new_tsj.append(tsj[x])
            # print ("LEN AFTER REMOVING NANS", len(new_tsi))
            #        print "Not a nan"
            # remove parts of the time series which are identical
            new_tsi = np.array(new_tsi)
            new_tsj = np.array(new_tsj)
            dts = new_tsi - new_tsj
            l = 0
            for x in range(len(dts)):
                if dts[x] == 0:
                    l += 1
            if l > 0:
                new_tsi, new_tsj = new_tsi[:l], new_tsj[:l]
            #print ("LEN OF UNIQUE SECTION: ", len(new_tsi))
            # take the log of them to test the correlation effect
            #log_tsi = np.log(new_tsi)
            #log_tsj = np.log(new_tsj)
            #cc[k] = np.corrcoef(new_tsi, new_tsj)[0, 1]
            #cc[k] = np.corrcoef(log_tsi, log_tsj)[0, 1]
            cc[k] = find_difference_between_arrays(new_tsi, new_tsj)
            k += 1

    #print np.isnan(cc)
    print cc
    # distances
    dd = np.arccos(cc)
    #print dd
    #print len(dd)
    # do agglomerative clustering by stepwise pair matching
    # based on angle between scalar products of time series
    ln = linkage(dd, method=method)

    # make a plot of the distance vs number of clusters. Use this to determine
    # the threshold
    thr = PlotDistanceVsNClusters(ln)

    # define threshold for cluster determination
    #thr = np.arccos(min_corr)
    #thr = 1.2

    # compute cluster indices
    cl = fcluster(ln, thr, criterion = 'distance')
    print("I've finished! I found {} clusters for you :)".format(cl.max()))
    print len(cl), n

    # assign the cluster id to the dataframe
    for i,id in enumerate(sources):
        df.loc[df.id==id, 'cluster_id'] = cl[i]

    # set colour palette: 8 class Set 1 from http://colorbrewer2.org
    N_colors = 8
    colors = LSDP.colours.list_of_hex_colours(N_colors, 'Dark2')[:cl.max()]
    threshold_color = '#377eb8'

    # now find the order of the cluster ids and assign the colours accordingly
    max_df = df.sort_values('distance_from_outlet', ascending=False).drop_duplicates(['cluster_id'])
    lengths = max_df['distance_from_outlet'].tolist()
    clusters = max_df['cluster_id'].tolist()
    print lengths
    print clusters
    sorted_clusters = [x for _,x in sorted(zip(lengths,clusters))]
    sorted_colors = []
    for i, c in enumerate(sorted_clusters):
        sorted_colors.append(colors[i])
        df.loc[df.cluster_id==c, 'colour'] = colors[i]

    set_link_color_palette(sorted_colors)

    plt.title('Hierarchical Clustering Dendrogram')
    #plt.xlabel('sample index')
    plt.ylabel('distance')
    R = dendrogram(ln, color_threshold=thr+0.00001, above_threshold_color=threshold_color, no_labels=True)

    plt.axhline(y = thr, color = 'r', ls = '--')
    plt.savefig(DataDirectory+fname_prefix+"_upstream_dendrogram.png", dpi=300)
    plt.clf()

    df.to_csv(DataDirectory+args.fname_prefix+'_profiles_upstream_clustered.csv', index=False)
    return df

def PlotDistanceVsNClusters(ln):
    """
    Make a plot of the distance between each cluster compared to the
    number of clusters. Maybe use this to determine where to put the distance
    threshold

    Args:
        ln: linkage matrix from clustering

    Author: FJC
    """
    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    #print ln
    # each iteration merges one cluster. so we start with n_clusters = n samples,
    # and then it reduces by one each time.
    clusters = []
    n_clusters = len(ln)+1
    for l in ln:
        # the distance is the 3rd column.
        this_dist = l[2]
        ax.scatter(n_clusters, this_dist, c='k')
        clusters.append(n_clusters)
        n_clusters -= 1

    # find the difference in the distances between each point in the linkage array
    dist = ln[:,2]
    deltas = [j-i for i, j in zip(dist[:-1], dist[1:])]
    # get the argmax of the difference
    i = np.argmax(deltas)
    n_clusters = clusters[i+1]
    # now find the distance threshold corresponding to this
    thr = dist[i]
    print ('The optimum distance threshold is '+str(thr))
    # now save the figure
    ax.set_xlabel('Number of clusters')
    ax.set_ylabel('Distance between clusters')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=2))
    plt.savefig(DataDirectory+fname_prefix+'_clusters_dist.png', dpi=300)
    plt.clf()

    return thr


def CalculateSlope(df, slope_window_size):
    """
    This function takes in a dataframe with elevation and distance
    from source data, and adds a column with the slope of each point
    fitted over a certain window size.

    Args:
        df: dataframe
        slope_window_size (int): total number of points used to calculate
        slope (INCLUDES the node of interest)

    Author: FJC
    """

    # make a plot of the gradient vs. distance from source aligned to the outlet
    fig_width_inches = 4.92126

    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # get all unique source IDs from the dataframe
    ids = df['id'].unique()
    for id in ids:
        dist = df['distance_from_outlet'][df['id'] == id]
        elev = df['elevation'][df['id'] == id]
        slopes = np.empty(len(dist))

        pts_array = np.column_stack((dist,elev))

        slicer = (slope_window_size - 1)/2

        for index, x in enumerate(pts_array):
            start_idx = index-slicer
            if start_idx < 0:
                start_idx=0
            end_idx = index+slicer+1
            if end_idx > len(pts_array):
                end_idx = len(pts_array)
            # find the rows above and below relating to the window size. We use whatever nodes
            # are available to not waste the data.
            this_slice = pts_array[start_idx:end_idx]
            # now regress this slice
            x = this_slice[:,0]
            y = this_slice[:,1]
            #print x, y
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            slopes[index] = abs(slope)
            #print slope

        df.loc[df.id==id, 'slope'] = slopes
        ax.plot(dist, slopes, lw=1)

    # now save the figure
    ax.set_xlabel('Distance from outlet (m)')
    ax.set_ylabel('Gradient')
    plt.savefig(DataDirectory+fname_prefix+'_profiles_upstream.png', dpi=300)
    plt.clf()

    print("Got the slope over a window radius of {} m").format(slope_window_size)

    return df

def RemoveProfilesShorterThanThresholdLength(df, profile_len=100):
    """
    Remove any profiles that are shorter than a threshold length.
    I can't believe I worked out how to do that in one line of code...!
    """
    # find the maximum distance for each id and remove any less than the profile length
    df = df[df.groupby('id')['distance_from_outlet'].transform('max') >= profile_len]
    return df

def RemoveProfilesWithShortUniqueSection(df, threshold_len=4):
    """
    Remove any profiles which only have a unique section that is
    shorter than a threshold length.
    """
    print ("Removing short unique profiles, the threshold length is: "+str(threshold_len))
    all_nodes = df['node'].tolist()
    sources = df['id'].unique()
    unique_sources = []
    for src in sources:
        this_df = df[df['id'] == src]
        these_nodes = this_df['node'].tolist()
        # check if each node is in the all nodes dataframe twice.
        # if it is then it is a duplicate.
        unique_nodes = 0
        for node in these_nodes:
            count = all_nodes.count(node)
            #print count
            if count < 2:
                unique_nodes += 1
        #print unique_nodes
        if unique_nodes > threshold_len:
            unique_sources.append(src)

    print ('Number of new sources: '+str(len(unique_sources)), 'number of old sources: '+str(len(sources)))
    df_new = df[df['id'].isin(unique_sources)]
    return df_new

def RemoveNonUniqueProfiles(df):
    """
    From the regularly spaced distance dataframe, remove any profiles
    which are non-unique (e.g., they have the same nodes as another
    profile). We check the first 5 upstream nodes of each profile and
    remove any which have a duplicate node within this.
    """
    # get a list of the sources
    sources = df['id'].unique()

    duplicate_sources = []
    # find the channel head of each source ID and check if this node is repeated
    # anywhere else in the dataframe
    for src in sources:
        if not src in duplicate_sources:
            this_df = df[df['id'] == src]
            rest_df = df[df['id'] != src]
        #     us_nodes = this_df['node'][-30:].tolist()
        #     for j in sources:
        #         # check the last 5 nodes of the rest of the dataframe
        #         rest_us_nodes = rest_df[rest_df['id'] == j]['node'][-30:].tolist()
        #         duplicates = filter(set(us_nodes).__contains__, rest_us_nodes)
        #         if duplicates:
        #             duplicate_sources.append(j)

            ch_node = this_df['node'].iloc[-1]
            mask = rest_df['node'].isin([ch_node]).any()
            if mask:
                duplicate_sources.append(src)

    df_new = df[~df['id'].isin(duplicate_sources)]
    return df_new
#---------------------------------------------------------------------#
# PLOTTING FUNCTIONS
#---------------------------------------------------------------------#
def PlotProfilesByCluster():
    """
    Function to make plots of the river profiles in each cluster

    Args:
        slope_window_size: window size in metres over which to calculate slope
        profile_len (int): number of data points in each profile (= distance from source)
        step (int): the spacing in metres between the data points that you want, default = 1m
        min_corr (float): minimum correlation threshold for clustering
        method (str): method for performing the clustering

    Author: FJC
    """

    # find the maximum distance from outlet in each cluster and use this to sort
    # the data
    # idx = cluster_df.groupby(['cluster_id'])['distance_from_outlet'].transform(max) == df['distance_from_outlet']
    # lengths = cluster_df[idx]
    # lengths = lengths.sort_values(by=['distance_from_outlet'])
    # print lengths
    # # find the unique clusters for plotting
    # clusters = lengths['cluster_id'].tolist()
    cluster_df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv')
    clusters = cluster_df['cluster_id'].unique()
    #clusters.sort()


    #counter = 0
    for cl in clusters:
        # set up a figure
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
        ax = fig.add_subplot(gs[5:100,10:95])

        this_df = cluster_df[cluster_df['cluster_id'] == cl]
        cl = int(this_df.iloc[0]['cluster_id'])
        this_colour = str(this_df.colour.unique()[0])
        sources = this_df['id'].unique()
        if (len(sources) > 1):
            for idx, src in enumerate(sources):
                src_df = this_df[this_df['id'] == src]
                src_df = src_df[src_df['slope'] != np.nan]
                ax.plot(src_df['distance_from_outlet'].as_matrix(), src_df['elevation'].as_matrix(), lw=1, color=this_colour)
                # save the colour to the cluster dataframe for later plots
                #cluster_df.loc[cluster_df.cluster_id==cl, 'colour'] = colors[counter]
            #counter +=1
        # else:
        #     ax.plot(this_df['distance_from_outlet'].as_matrix(), this_df['elevation'].as_matrix(), lw=1, color=threshold_color)
            # save the colour to the cluster dataframe for later plots
            #cluster_df.loc[cluster_df.cluster_id==cl, 'colour'] = threshold_color

        ax.set_xlabel('Distance from outlet (m)')
        ax.set_ylabel('Elevation (m)')
        ax.set_title('Cluster {}'.format(int(cl)))

        plt.savefig(DataDirectory+fname_prefix+('_profiles_upstream_clustered_{}.png').format(int(cl)), dpi=300)
        plt.clf()

    # write the clustered dataframe to csv
    #cluster_df.to_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv', index=False)

    return cluster_df

def MakeHillshadePlotClusters(drape_lithology=False):
    """
    Make a shaded relief plot of the raster with the channels coloured by the cluster
    value. Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

    Args:
        drape_lithology (bool): true to add in the raster with the lithology information. At the
        moment only works for the model runs.

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
    cluster_df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv')


    # set figure sizes based on format
    fig_width_inches = 4.92126

    # some raster names
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext

    # create the map figure
    if drape_lithology:
        MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM",colourbar_location='bottom')
        LithoName = 'spatial_K_KRaster'+raster_ext
        MF.add_drape_image(LithoName, DataDirectory,alpha=0.5,colourmap="gray", show_colourbar = True, colorbarlabel='K', discrete_cmap=True, n_colours=2,cbar_type=str)
    else:
        MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM")
    clusters = cluster_df.cluster_id.unique()
    for cl in clusters:
        # plot the whole channel network in black
        #ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
        #MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='0.8',manual_size=5, zorder=1)
        # plot the clustered profiles in the correct colour
        this_df = cluster_df[cluster_df.cluster_id == cl]
        this_colour = str(this_df.colour.unique()[0])
        ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
        MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=5)

    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = DataDirectory+fname_prefix+'_hs_clusters.png', FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False) # Save the figure

def PlotMedianProfiles():
    """
    Make a summary plot showing the median profile for each cluster, both in
    gradient-distance and elevation-distance space.

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv')

    # find out some info
    clusters = df.cluster_id.unique()
    clusters.sort()
    sources = df.id.unique()
    dist_array = df[df.id == sources[0]].reg_dist.as_matrix()

    # set up a figure
    fig,ax = plt.subplots(nrows=len(clusters),ncols=1, figsize=(5,6), sharex=True, sharey=True)
    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    # for each cluster, get the mean gradient for each regular distance
    for i, cl in enumerate(clusters):

        cluster_df = df[df.cluster_id == cl]
        median_gradients = np.asarray([cluster_df[cluster_df.reg_dist == x].slope.median() for x in dist_array])
        lower_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].slope.quantile(0.25) for x in dist_array])
        upper_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].slope.quantile(0.75) for x in dist_array])
        # get the colour from the dataframe
        this_colour = str(cluster_df.colour.unique()[0])
        ax[i].plot(dist_array,median_gradients,color=this_colour, lw=1)
        ax[i].fill_between(dist_array, lower_quantile, upper_quantile, facecolor=this_colour, alpha=0.2)
        ax[i].text(0.9, 0.8,'Cluster {}'.format(int(cl)),horizontalalignment='center',verticalalignment='center',transform = ax[i].transAxes,fontsize=8)

    # set axis labels
    plt.xlabel('Distance from outlet (m)')
    plt.ylabel('Gradient')

    # save and clear the figure
    plt.savefig(DataDirectory+fname_prefix+('_profiles_median.png'), dpi=300)
    plt.clf()
    plt.cla()
    plt.close()

    # set up a figure
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # for each cluster, get the mean gradient for each regular distance
    for cl in clusters:
        cluster_df = df[df.cluster_id == cl]
        median_elevs = np.asarray([cluster_df[cluster_df.reg_dist == x].elevation.median() for x in dist_array])
        lower_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].elevation.quantile(0.25) for x in dist_array])
        upper_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].elevation.quantile(0.75) for x in dist_array])
        # get the colour from the dataframe
        this_colour = str(cluster_df.colour.unique()[0])
        ax.plot(dist_array,median_elevs,color=this_colour, lw=1)
        ax.fill_between(dist_array, lower_quantile, upper_quantile, facecolor=this_colour, alpha=0.2)

    ax.set_xlabel('Distance from outlet (m)')
    ax.set_ylabel('Elevation (m)')

    plt.savefig(DataDirectory+fname_prefix+('_profiles_median_elev.png'), dpi=300)
    plt.clf()

def PlotSlopeArea():
    """
    Make a summary plot showing the S-A plot for each cluster.

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv')

    # find out some info
    clusters = df.cluster_id.unique()
    clusters.sort()
    sources = df.id.unique()

    # set up a figure
    fig,ax = plt.subplots(nrows=len(clusters),ncols=1, figsize=(5,8), sharex=False, sharey=False)
    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    # for each cluster, get the mean gradient for each regular distance
    for i, cl in enumerate(clusters):

        cluster_df = df[df.cluster_id == cl]

        # get the colour from the dataframe
        this_colour = str(cluster_df.colour.unique()[0])
        ax[i].scatter(cluster_df['drainage_area'], cluster_df['slope'], color=this_colour, s=1)
        ax[i].text(0.9, 0.9,'Cluster {}'.format(int(cl)),horizontalalignment='center',verticalalignment='center',transform = ax[i].transAxes,fontsize=8)
        ax[i].set_xscale('log')
        ax[i].set_yscale('log')
        ax[i].set_ylim(0.0001, 1)

    # set axis labels
    plt.xlabel('Drainage area (m$^2$)')
    plt.ylabel('Gradient', labelpad=15)
    plt.subplots_adjust(left=0.15, hspace=0.5)

    # save and clear the figure
    plt.savefig(DataDirectory+fname_prefix+('_SA_median.png'), dpi=300)
    plt.clf()
    plt.cla()
    plt.close()

def PlotUniqueStreamsWithLength(step=2, slope_window_size=25):
    """
    Function to make a plot of the number of unique channels you get (at least one non-overlapping
    node with different profile lengths that you analyse).

    Author: FJC
    """
    # read in the original csv
    df = pd.read_csv(DataDirectory+args.fname_prefix+'_all_tribs.csv')

    max_length = df['distance_from_outlet'].max()
    len_step=100
    profile_lengths = np.arange(len_step, (max_length/100).astype(int)*100, step=len_step)

    # calculate the slope
    df = CalculateSlope(df, slope_window_size)

    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # for each length, remove any sources that are below the desired length.
    for len in profile_lengths:
        thinned_df = RemoveProfilesShorterThanThresholdLength(df, len)
        reg_df, __ = ProfilesRegularDistance(thinned_df, profile_len = len, step=step, slope_window_size=slope_window_size)
        reg_df = RemoveNonUniqueProfiles(reg_df)
        n_sources = reg_df['id'].unique().size
        print n_sources
        ax.scatter(len, n_sources, s = 25, c='w', edgecolors='k')

    ax.set_xlabel('Profile length (m)')
    ax.set_ylabel('Number of unique channels')

    plt.savefig(DataDirectory+fname_prefix+'_n_channels_with_length.png', dpi=300)
    plt.clf()

def PlotLongitudinalProfiles():
    """
    Just make a simple plot of the river long profiles
    """
    df = pd.read_csv(DataDirectory+args.fname_prefix+'_profiles_upstream_clustered.csv')

    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    sources = df['id'].unique()

    for src in sources:
        this_df = df[df['id'] == src]
        ax.plot(this_df['distance_from_outlet'], this_df['elevation'])

    ax.set_xlabel('Distance from outlet (m)')
    ax.set_ylabel('Elevation (m)')

    plt.savefig(DataDirectory+fname_prefix+'_long_profiles.png', dpi=300)
    plt.clf()

def MakeShadedSlopeMap():
    """
    Make a nice shaded slope image
    """
    # set figure sizes based on format
    fig_width_inches = 4.92126

    # some raster names
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    SlopeName = fname_prefix+'_slope'+raster_ext

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM",colourbar_location='right')
    MF.add_drape_image(BackgroundRasterName, DataDirectory,alpha=0.4,colourmap="terrain", show_colourbar = True, colorbarlabel='Elevation (m)', discrete_cmap=False)

    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = DataDirectory+fname_prefix+'_hs_slope.png', FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False) # Save the figure

def PlotTrunkChannel():
    """
    Make a simple plot of the longest channel. This is mostly to use for the model runs.
    """
    df = pd.read_csv(DataDirectory+args.fname_prefix+'_profiles_upstream_clustered.csv')

    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    trunk_src = df.loc[df['distance_from_outlet'].idxmax()]['id']

    this_df = df[df['id'] == trunk_src]
    ax.plot(this_df['distance_from_outlet'], this_df['elevation'], c='k')

    ax.set_xlabel('Distance from outlet (m)')
    ax.set_ylabel('Elevation (m)')
    #ax.set_xlim(0,2500)
    #ax.set_ylim(0,35)

    plt.savefig(DataDirectory+fname_prefix+'_trunk_profile.png', dpi=300)
    plt.clf()

if __name__ == '__main__':

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()

    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the m/n analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error (unless you're running the parallel plotting).")

    # The options for clustering
    parser.add_argument("-len", "--profile_len", type=int, help="The minimum length over which a profile has to be unique (i.e. does not share nodes with any other profile)", default=10)
    parser.add_argument("-sw", "--slope_window", type=int, help="The window size for calculating the slope based on a regression through an equal number of nodes upstream and downstream of the node of interest. This is the total number of nodes that are used for calculating the slope. For example, a slope window of 25 would fit a regression through 12 nodes upstream and downstream of the node, plus the node itself. The default is 25 nodes.", default=25)
    parser.add_argument("-m", "--method", type=str, help="The method for clustering, see the scipy linkage docs for more information. The default is 'ward'.", default='ward')
    parser.add_argument("-c", "--min_corr", type=float, help="The minimum correlation for defining the clusters. Use a smaller number to get less clusters, and a bigger number to get more clusters (from 0 = no correlation, to 1 = perfect correlation). The default is 0.5. DEPRECATED - now we calculate the threshold statistically.", default=0.5)
    parser.add_argument("-step", "-step", type=int, help="The regular spacing in metres that you want the profiles to have for the clustering. This should be greater than sqrt(2* DataRes^2).  The default is 2 m which is appropriate for grids with a resolution of 1 m.", default = 2)

    args = parser.parse_args()

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()
    else:
        fname_prefix = args.fname_prefix

    # get the base directory
    if args.base_directory:
        DataDirectory = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not DataDirectory.endswith("/"):
            print("You forgot the '/' at the end of the directory, appending...")
            DataDirectory = DataDirectory+"/"
    else:
        print("WARNING! You haven't supplied the data directory. I'm using the current working directory.")
        DataDirectory = os.getcwd()

    # print the arguments that you used to an output file for reproducibility
    with open(DataDirectory+args.fname_prefix+'_report.csv', 'w') as output:
        for arg in vars(args):
            output.write(str(arg)+','+str(getattr(args, arg))+'\n')
        output.close()

    # check if the slopes file exists
    slope_file = DataDirectory+args.fname_prefix+'_slopes.csv'
    if os.path.isfile(slope_file):
        df = pd.read_csv(slope_file)
    else:
        # read in the original csv
        df = pd.read_csv(DataDirectory+args.fname_prefix+'_all_tribs.csv')
        #git stdf = RemoveProfilesWithShortUniqueSection(df, args.profile_len)
        # calculate the slope
        df = CalculateSlope(df, args.slope_window)
        df.to_csv(DataDirectory+args.fname_prefix+'_slopes.csv', index=False)

    cluster_df = ClusterProfilesDrainageArea(df, profile_len = args.profile_len, step=args.step)
    # regular_df = ProfilesRegDistVaryingLength(df, profile_len=args.profile_len, step=args.step, slope_window_size=args.slope_window)
    #
    # # cluster the profiles
    # regular_df = pd.read_csv(DataDirectory+args.fname_prefix+'_profiles_upstream_reg_dist_var_length.csv')
    # cluster_df = ClusterProfilesVaryingLength(regular_df, args.min_corr,method=args.method)
    PlotProfilesByCluster()
    # #
    # #PlotMedianProfiles()
    MakeHillshadePlotClusters()
    PlotSlopeArea()
    PlotTrunkChannel()
    # # PlotLongitudinalProfiles()
    #MakeShadedSlopeMap()

    print('Enjoy your clusters, pal')

    #PlotUniqueStreamsWithLength()
