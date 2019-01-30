import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
from scipy import stats
import math
import plotting as pl

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
    #print norm
    #print np.sqrt(n)
    diff = 1-(norm/np.sqrt(n))
    #diff = 1-norm
    #print diff
    return diff

def AverageEuclidianDifference(x, y):
    """
    Find the average Euclidian difference between two arrays x and y:
    d = sqrt(sum(x-y)^2)/n
    Liao, 2005,  doi:10.1016/j.patcog.2005.01.025
    """
    n = len(x)
    d = (np.sqrt(np.sum((x - y)*(x - y))))/n
    return d

#def MinimiseLag(x, y, s=200):
    """
    Take in two arrays x and y, and do a shifting to minimise the average
    euclidian distance between them

    Args:
        x: first array
        y: second array
        s: maximum shift, default = 200 px
    FJC
    """



def ProfilesRegularDistance(DataDirectory, fname_prefix, df, profile_len = 1000, step=2, slope_window_size=25):
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
            slopes = this_df['slope'].values[::-1] #need to reverse these as the distances need to be sorted
            distances = this_df['distance_from_outlet'].values[::-1]
            if (distances.max() >= profile_len):
                profiles.append((distances, slopes))
                #thinned_df = thinned_df.append(this_df)
                final_sources.append(source)

    # now create the 2d array to store the data
    n_profiles = len(profiles)
    data = np.empty((n_profiles, len(reg_dist)))

    # create a new dataframe for storing the data about the selected profiles
    thinned_df = pd.DataFrame()

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

    # write the thinned_df to output in case we want to reload
    thinned_df.to_csv(DataDirectory+fname_prefix+'_profiles_upstream_reg_dist.csv',index=False)
    #
    return thinned_df, data

def ProfilesRegDistVaryingLength(DataDirectory, fname_prefix, df, profile_len=4,step=2, slope_window_size=25):
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
    print (source_ids)
    rows_list = []
    for i, source in enumerate(source_ids):
        this_df = df[df['id'] == source]
        # create new array of regularly spaced differences
        reg_dist = np.arange(step, int(this_df.distance_from_outlet.max())+step, step)

        if not this_df.empty:
            df_array = this_df.values[::-1]
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

def GetProfilesByStreamOrder(DataDirectory, fname_prefix, df,step=2,slope_window_size=25,stream_order=1):
    """
    Take the dataframe and return only the profiles of a certain stream order,
    regularly sampled.
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

    # only keep the nodes belonging to the corresponding stream order
    so_df = df[df['stream_order'] == stream_order]
    longest_df = so_df

    # if not first order streams then find the longest channel
    if stream_order != 1:
        # get a list of the first highest order node for checking
        check_nodes = so_df.groupby(by='id').apply(lambda x: x.iloc[-1]['node']).tolist()

        # now find the longest channel in the df that has this node in it
        so_df = df[df['stream_order'] <= stream_order]
        new_df = so_df.groupby(by='id')
        longest_ids = []
        # loop through each higher order node and find all the channel profiles that contain them
        for node in check_nodes:
            these_ids = new_df.apply(lambda x: x[x['node'] == node])['id'].tolist()   # get the ids that have this node in it
            this_df = so_df[so_df['id'].isin(these_ids)]  # filter the SO dataframe to only these IDs
            # now find the ID of the longest one
            longest_ids.append(this_df.loc[this_df['distance_from_outlet'].idxmax()]['id'])
        #print(longest_ids)
        longest_df = so_df[so_df['id'].isin(longest_ids)]

    #print so_df

    # loop through the dataframe and store the data for each profile as an array of
    # slopes and distances
    source_ids = longest_df['id'].unique()
    #print (source_ids)
    rows_list = []
    for i, source in enumerate(source_ids):
        this_df = longest_df[longest_df['id'] == source]
        # get the distances from the channel head
        distances = [this_df.distance_from_outlet.max()-x for x in this_df.distance_from_outlet]
    #    print distances
        # create new array of regularly spaced differences
        reg_dist = np.arange(step, int(np.max(distances)+step), step)

        if not this_df.empty:
            df_array = this_df.values
            slopes = this_df['slope'].values
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

    # write the thinned_df to output in case we want to reload
    thinned_df.to_csv(DataDirectory+fname_prefix+'_profiles_SO{}.csv'.format(stream_order), index=False)

    # now save the figure
    ax.set_xlabel('Distance from source (m)')
    ax.set_ylabel('Gradient')
    plt.savefig(DataDirectory+fname_prefix+'_profiles_SO{}.png'.format(stream_order), dpi=300)
    plt.clf()

    return thinned_df

def ClusterProfiles(DataDirectory, fname_prefix, df, profile_len=100, step=2, min_corr=0.5, method='complete'):
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

    set_link_color_palette(colors)

    source_ids = df['id'].unique()

    plt.title('Hierarchical Clustering Dendrogram')
    plt.ylabel('distance')
    R = dendrogram(ln, color_threshold=thr, above_threshold_color=threshold_color,no_labels=True)

    plt.axhline(y = thr, color = 'r', ls = '--')
    plt.savefig(DataDirectory+fname_prefix+"_upstream_dendrogram.png", dpi=300)
    plt.clf()

    for i,id in enumerate(source_ids):
        df.loc[df.id==id, 'cluster_id'] = cl[i]

    return df

def ClusterProfilesVaryingLength(DataDirectory, OutDirectory, fname_prefix, df, method='ward',stream_order=1,threshold_level=0):
    """
    Cluster the profiles based on gradient and distance from source. This works for profiles of varying length.
    Aggolmerative clustering, see here for more info:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage

    Args:
        df: pandas dataframe from the river profile csv.
        method (str): clustering method to use, see scipy docs. Can be 'single', 'complete', 'average',
        'weighted', 'centroid', 'median', or 'ward'. Default is 'ward'.
        stream_order: the stream order of the profiles, default = first order.
        threshold_level: the level at which to cut the dendrogram. Threshold level 0 is the default, and this is the level
        with the maximum distance between clusters. Can increase this to 1 to take the second max distance.

    Author: AR, FJC
    """
    print ("Now I'm going to do some hierarchical clustering...")
    np.set_printoptions(threshold='nan')

    #sort the dataframe based on max distance from outlet for each source id.

    # get the data from the dataframe into the right format for clustering
    sources = df['id'].unique()
    n = len(sources)
    data = []

    for i, src in enumerate(sources):
        this_df = df[df['id'] == src]
        data.append(this_df['slope'].values)

    # correlation coefficients
    cc = np.zeros(int(n * (n - 1) / 2))
    #cc = []
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            tsi = data[i]
            tsj = data[j]
            if len(tsi) > len(tsj):
                tsi = tsi[:len(tsj)]
            else:
                tsj = tsj[:len(tsi)]
            cc[k] = AverageEuclidianDifference(tsi, tsj)
            k += 1

    ln = linkage(cc, method=method)

    # make a plot of the distance vs number of clusters. Use this to determine
    # the threshold
    thr = PlotDistanceVsNClusters(DataDirectory, fname_prefix, ln, threshold_level)

    # compute cluster indices
    cl = fcluster(ln, thr, criterion = 'distance')
    print("I've finished! I found {} clusters for you :)".format(cl.max()))
    #print([int(c) for c in cl])

    # assign the cluster id to the dataframe
    for i,src in enumerate(sources):
        df.loc[df.id==src, 'cluster_id'] = cl[i]

    # set colour palette: 8 class Set 1 from http://colorbrewer2.org
    N_colors = 8
    #colors = pl.list_of_hex_colours(N_colors, 'Set1')[:cl.max()]
    colors = ['#e41a1c', '#377eb8', '#000000']
    threshold_color = '#A9A9A9'
    clusters = df['cluster_id'].unique()

    # max_df = df.sort_values('distance_from_outlet', ascending=False).drop_duplicates(['cluster_id'])
    # lengths = max_df['distance_from_outlet'].tolist()
    # print lengths
    # print clusters
    sorted_clusters = sorted(clusters)
    # sorted_colors = [x for _,x in sorted(zip(lengths,colors))]
    for i, c in enumerate(sorted_clusters):
        df.loc[df.cluster_id==c, 'colour'] = colors[i]

    set_link_color_palette(colors)

    source_ids = df['id'].unique()

    #plt.title('Hierarchical Clustering Dendrogram')
    plt.ylabel('Dissimilarity ($d_R$)', fontsize=14)
    R = dendrogram(ln, color_threshold=thr+0.00001, above_threshold_color=threshold_color,no_labels=True)

    plt.axhline(y = thr, color = 'r', ls = '--')
    plt.savefig(OutDirectory+fname_prefix+"_dendrogram_SO{}.png".format(stream_order), dpi=300)
    plt.clf()

    df.to_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order), index=False)

    return df

def ClusterProfilesDrainageArea(DataDirectory, fname_prefix, df, profile_len=100, step=2, method='ward'):
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

    all_areas = df['drainage_area'].values
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

        #df_array = this_df.values[::-1]
        slopes = this_df['slope'].values
        areas = this_df['drainage_area'].values
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
            new_tsi, new_tsj = new_tsi[:l], new_tsj[:l]
            #print ("LEN OF UNIQUE SECTION: ", len(new_tsi))
            # take the log of them to test the correlation effect
            log_tsi = np.log(new_tsi)
            log_tsj = np.log(new_tsj)
            #cc[k] = np.corrcoef(new_tsi, new_tsj)[0, 1]
            #cc[k] = np.corrcoef(log_tsi, log_tsj)[0, 1]
            cc[k] = find_difference_between_arrays(new_tsi, new_tsj)
            k += 1

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
    print (len(cl), n)

    # assign the cluster id to the dataframe
    for i,id in enumerate(sources):
        df.loc[df.id==id, 'cluster_id'] = cl[i]

    # set colour palette: 8 class Set 1 from http://colorbrewer2.org
    N_colors = 8
    colors = list_of_hex_colours(N_colors, 'Dark2')[:cl.max()]
    threshold_color = '#377eb8'

    # now find the order of the cluster ids and assign the colours accordingly
    max_df = df.sort_values('distance_from_outlet', ascending=False).drop_duplicates(['cluster_id'])
    lengths = max_df['distance_from_outlet'].tolist()
    clusters = max_df['cluster_id'].tolist()
    print (lengths)
    print (clusters)
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

def PlotDistanceVsNClusters(DataDirectory, fname_prefix, ln, threshold_level=0):
    """
    Make a plot of the distance between each cluster compared to the
    number of clusters. Maybe use this to determine where to put the distance
    threshold

    Args:
        ln: linkage matrix from clustering
        threshold_level: which level to return the threshold at. 0 = max distance between clusters.

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
    n_clusters = clusters[i-threshold_level]
    # now find the distance threshold corresponding to this
    thr = dist[i-threshold_level]
    print ('The optimum distance threshold is '+str(thr))
    # now save the figure
    plt.xlabel('Number of clusters')
    plt.ylabel('Distance between clusters')
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(base=2))
    plt.savefig(DataDirectory+fname_prefix+'_clusters_dist.png', dpi=300)
    plt.clf()

    return thr


def CalculateSlope(DataDirectory, fname_prefix, df, slope_window_size):
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
            this_slice = pts_array[int(start_idx):int(end_idx)]
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

    print("Got the slope over a window radius of {} m".format(slope_window_size))

    return df

def RemoveProfilesShorterThanThresholdLength(df, profile_len=5):
    """
    Remove any profiles that are shorter than a threshold length.
    (number of nodes)
    """
    # find the maximum distance for each id and remove any less than the profile length
    df = df.loc[df.groupby('id').filter(lambda x: len(x) >= profile_len).index]
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
    profile).
    """
    duplicate_sources = []
    # get a list of the sources
    sources = df['id'].unique()
    n = len(sources)
    for i in range(n):
        for j in range(i+1, n):
            df1 = df[df.id == sources[i]].node.tolist()
            df2 = df[df.id == sources[j]].node.tolist()
            if not set(df1).isdisjoint(df2):
                duplicate_sources.append(sources[j])

    df_new = df[~df['id'].isin(duplicate_sources)]
    return df_new
#---------------------

def PrintJunctionsToCSV(DataDirectory, fname_prefix, stream_order=1):
    """
    print a list of the downstream junction of the channels in each cluster
    for extracting the catchments using lsdtopotools
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))
    clusters = df.cluster_id.unique()

    for cl in clusters:
        this_df = df[df.cluster_id == cl]
        sources = this_df.id.unique()
        ds_nodes = this_df.loc[this_df.groupby('id')['reg_dist'].idxmax()]
        ds_nodes.to_csv(DataDirectory+fname_prefix+'_junctions_SO{}_CL{}.csv'.format(stream_order, int(cl)), index=False)
