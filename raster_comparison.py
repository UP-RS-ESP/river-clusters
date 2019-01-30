# Scripts to compare the results of the clustering to different rasters
# FJC 12/09/18

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from osgeo import ogr
from shapely.geometry import shape, Polygon
from descartes.patch import PolygonPatch
import matplotlib.cm as cm
import LSDPlottingTools as LSDP
from LSDPlottingTools import LSDMap_GDALIO as IO
from LSDPlottingTools import LSDMap_PointTools as PT
import os

def BoxPlotByCluster(DataDirectory, OutDirectory, fname_prefix,  raster_name, stream_order=1):
    """
    Make a boxplot of the results of the clustering compared to the raster specified
    by raster_name
    """
    #df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))

    # read in the raster
    raster_ext = '.bil'
    this_raster = IO.ReadRasterArrayBlocks(DataDirectory+raster_name)
    EPSG_string = IO.GetUTMEPSG(DataDirectory+raster_name)
    NDV, xsize, ysize, GeoT, Projection, DataType = IO.GetGeoInfo(DataDirectory+raster_name)
    CellSize,XMin,XMax,YMin,YMax = IO.GetUTMMaxMin(DataDirectory+raster_name)

    pts = PT.LSDMap_PointData(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order),data_type ='csv')
    easting, northing = pts.GetUTMEastingNorthing(EPSG_string=EPSG_string)
    cluster_id = pts.QueryData('cluster_id', PANDEX=True)
    clusters = list(set(cluster_id))

    # dict for the data
    data = {k: [] for k in clusters}
    for x, (i, j) in enumerate(zip(northing, easting)):
    # convert to rows and cols
        X_coordinate_shifted_origin = j - XMin;
        Y_coordinate_shifted_origin = i - YMin;

        col_point = int(X_coordinate_shifted_origin/CellSize);
        row_point = (ysize - 1) - int(round(Y_coordinate_shifted_origin/CellSize));
        # check for data at this cell
        this_value = this_raster[row_point][col_point]
        if not np.isnan(this_value):
            if this_value < 10:
                # get the cluster id
                data[cluster_id[x]].append(this_value)

    print(data)

    # now make a boxplot
    labels, these_data = [*zip(*data.items())]  # 'transpose' items to parallel key, value lists

    plt.boxplot(these_data)
    plt.xticks(range(1, len(labels) + 1), labels)
    plt.show()

def BoxPlotGradient(DataDirectory, OutDirectory, fname_prefix, stream_order=1):
    """
    Make a boxplot showing the channel gradient stats for each cluster
    """

    # read the csv and get some info
    df = pd.read_csv(OutDirectory+fname_prefix+"_profiles_clustered_SO{}.csv".format(stream_order))
    colors = df['colour'].unique()

    # set props for fliers
    flierprops = dict(marker='o', markerfacecolor='none', markersize=1,
                  linestyle='none', markeredgecolor='k')

    # make the boxplot and return the dict with the boxplot properties
    bp_dict = df.boxplot(column=['slope'], by=['cluster_id'], return_type='both', patch_artist= True,flierprops=flierprops, figsize=(5,5))
    # make the median lines black
    #[[item.set_color('k') for item in bp_dict[key]['medians']] for key in bp_dict.keys()]

    # change the colours based on the cluster ID
    for row_key, (ax,row) in bp_dict.iteritems():
        plt.xlabel('')
        ax.set_title('')
        j=-1 #stupid thing because there are double the number of caps and whiskers compared to boxes
        for i,cp in enumerate(row['caps']):
            if i%2==0:
                j+=1
            cp.set(color=colors[j])
        j=-1
        for i,wh in enumerate(row['whiskers']):
            if i%2==0:
                j+=1
            wh.set_color(colors[j])
        for i,box in enumerate(row['boxes']):
            box.set_facecolor(colors[i])
            box.set_alpha(0.7)
            box.set_edgecolor(colors[i])
        for i,med in enumerate(row['medians']):
            med.set(color=colors[i])
        for i,pt in enumerate(row['fliers']):
            pt.set_markeredgecolor(colors[i])


    ax.grid(color='0.8', linestyle='--', which='major', zorder=1)
    ax.set_xticklabels(['Canada shale', 'Breccia/volcaniclastics'], fontsize=12)
    #print(boxplot)
    plt.suptitle('')
    ax.set_ylabel('Gradient (m/m)', fontsize=14)
    plt.subplots_adjust(left=0.2)
    plt.savefig(OutDirectory+fname_prefix+'_boxplot.png', dpi=300)
    plt.clf()

def GetLithologyPercentages(DataDirectory, OutDirectory, fname_prefix, raster_name, stream_order=1):
    """
    Get the percentage of the nodes in each cluster that drain each lithology
    """
    from collections import Counter
    # read in the raster
    raster_ext = '.bil'
    this_raster = IO.ReadRasterArrayBlocks(DataDirectory+raster_name)
    EPSG_string = IO.GetUTMEPSG(DataDirectory+raster_name)
    NDV, xsize, ysize, GeoT, Projection, DataType = IO.GetGeoInfo(DataDirectory+raster_name)
    CellSize,XMin,XMax,YMin,YMax = IO.GetUTMMaxMin(DataDirectory+raster_name)

    pts = PT.LSDMap_PointData(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order),data_type ='csv')
    easting, northing = pts.GetUTMEastingNorthing(EPSG_string=EPSG_string)
    cluster_id = pts.QueryData('cluster_id', PANDEX=True)
    clusters = list(set(cluster_id))

    # dict for the data
    data = {k: [] for k in clusters}
    for x, (i, j) in enumerate(zip(northing, easting)):
    # convert to rows and cols
        X_coordinate_shifted_origin = j - XMin;
        Y_coordinate_shifted_origin = i - YMin;

        col_point = int(X_coordinate_shifted_origin/CellSize);
        row_point = (ysize - 1) - int(round(Y_coordinate_shifted_origin/CellSize));
        # check for data at this cell
        this_value = this_raster[row_point][col_point]
        if not np.isnan(this_value):
            data[cluster_id[x]].append(this_value)

    # you have the values. now what percentage are each?
    for key, liths in data.items():
        c = Counter(liths)
        n_ndv = c[0.0]
        print(c)
        [print(x,": ",vals/len(liths) * 100) for x, vals in c.items()]

def ReadBasinPolygons(DataDirectory, OutDirectory, raster_name):
    """
    Read in the basin polygons
    """
    import shapefile

    #ax.set_aspect('equal')

    # check if the shapefile exists
    shpfile = raster_name+'.shp'
    if not os.path.isfile(OutDirectory+shpfile):
        print("Polygonising the basin raster...")
        # read in the raster
        raster_ext = '.bil'
        polygons = IO.PolygoniseRaster(OutDirectory, raster_name+raster_ext, raster_name)
        polygons = list(polygons.values())

    else:
        # read in the shapefile
        sf  = shapefile.Reader(OutDirectory+shpfile)
        # shape = file.GetLayer(0)

        polygons = []
        for s in list(sf.iterShapes()):
            nparts = len(s.parts) # total parts
            if nparts == 1:
               polygon = Polygon(s.points)
               polygons.append(polygon)

            else: # loop over parts of each shape, plot separately
              for ip in range(nparts): # loop over parts, plot separately
                  i0=s.parts[ip]
                  if ip < nparts-1:
                     i1 = s.parts[ip+1]-1
                  else:
                     i1 = len(s.points)

                  polygon = Polygon(s.points[i0:i1+1])
                  polygons.append(polygon)

    return polygons

def PlotBasinsWithHillshade(DataDirectory, OutDirectory, fname_prefix, stream_order=1):
    """
    Read in the basins and plot them over a hillshade coloured by their cluster ID
    """

    df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))
    clusters = df.cluster_id.unique()

    # make a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.1,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,5:100])

    # plot the raster
    hs_raster = IO.ReadRasterArrayBlocks(DataDirectory+fname_prefix+'_hs.bil')
    extent = LSDP.GetRasterExtent(DataDirectory+fname_prefix+'_hs.bil')
    plt.imshow(hs_raster, cmap=cm.gray, extent=extent)

    means = {}
    for i,cl in enumerate(clusters):
        this_df = df[df.cluster_id == cl]
        # get the polygons
        polygons = ReadBasinPolygons(DataDirectory, OutDirectory, fname_prefix+'_basins_SO{}_CL{}'.format(stream_order,int(cl)))
        for p in polygons:
            #print(list(p.exterior.coords))
            patch = PolygonPatch(p, facecolor=this_df.iloc[0]['colour'], alpha=1.0, zorder=2, lw=0.2)
            ax.add_patch(patch)
        #print(polygons)
        # for each

    plt.savefig(OutDirectory+'polygons.png', FigFormat='png', dpi=500)

def MakeBoxPlotsKsnLithology(DataDirectory, fname_prefix, raster_name, theta=0.45, label_list=[]):
    """
    Make boxplots of ksn compared to lithology raster. Lithology should have integer
    values for the different rock types (rock type with 0 is excluded). Pass in list of
    labels for the different units, which must be the same length as the number of lithology codes.
    If none is passed then just use the integer values for labelling.
    """
    from scipy import stats

    # read in the raster
    raster_ext = '.bil'
    this_raster = IO.ReadRasterArrayBlocks(DataDirectory+raster_name)
    #EPSG_string = IO.GetUTMEPSG(DataDirectory+raster_name)
    EPSG_string='epsg:32611'
    print(EPSG_string)
    NDV, xsize, ysize, GeoT, Projection, DataType = IO.GetGeoInfo(DataDirectory+raster_name)
    CellSize,XMin,XMax,YMin,YMax = IO.GetUTMMaxMin(DataDirectory+raster_name)

    pts = PT.LSDMap_PointData(DataDirectory+fname_prefix+'_ksn.csv',data_type ='csv')
    print(pts)
    easting, northing = pts.GetUTMEastingNorthing(EPSG_string=EPSG_string)
    ksn = pts.QueryData('ksn', PANDEX=True)
    #print(ksn)

    # get the unique values in the raster
    raster_values = np.unique(this_raster)
    raster_values = raster_values[1:]


    # dict for the data
    data = {k: [] for k in raster_values}

    for x, (i, j) in enumerate(zip(northing, easting)):
    # convert to rows and cols
        X_coordinate_shifted_origin = j - XMin;
        Y_coordinate_shifted_origin = i - YMin;

        col_point = int(X_coordinate_shifted_origin/CellSize);
        row_point = (ysize - 1) - int(round(Y_coordinate_shifted_origin/CellSize));
        # check for data at this cell
        this_raster_value = this_raster[row_point][col_point]
        if not np.isnan(this_raster_value) and this_raster_value != 0:
            data[this_raster_value].append(ksn[x])

    # set up a figure
    fig,ax = plt.subplots(nrows=1,ncols=1, figsize=(5,5), sharex=True, sharey=True)

    labels, dict = [*zip(*data.items())]  # 'transpose' items to parallel key, value lists
    print(label_list)
    if label_list:
        labels = label_list
    print(labels)
    box = plt.boxplot(dict, patch_artist=True)
    plt.xticks(range(1, len(labels) + 1), labels)
    plt.ylabel('$k_{sn}$', fontsize=14)

    # get the medians for plotting as an upper label
    medians = []
    print("========SOME KSN STATISTICS=========")
    for key, value in data.items():
        print("Key {}, median ksn = {}".format(key, np.median(value)))
        medians.append(np.median(value))
        print("Key {}, IQR = {}".format(key, stats.iqr(value)))
    print("========================================")
    pos = np.arange(len(labels)) + 1
    upperLabels = [str(np.round(s, 2)) for s in medians]

    # change the colours for each lithology
    colors=['#60609fff', '#fdbb7fff', '#935353ff', '#f07b72ff']
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.9)
        patch.set_edgecolor('k')
    for cap in box['caps']:
        cap.set(color='k')
    for wh in box['whiskers']:
        wh.set(color='k')
    for med in box['medians']:
        med.set(color='k')
    for flier, color in zip(box['fliers'], colors):
        flier.set_markeredgecolor(color)
        flier.set_markerfacecolor(color)
        flier.set_markersize(2)

    # for tick, label in zip(range(len(labels)), ax.get_xticklabels()):
    #     k = tick % 2
    #     ax.text(pos[tick], top - (top*0.05), upperLabels[tick],
    #              horizontalalignment='center', color=colors[k])

    ax.grid(color='0.8', linestyle='--', which='major', zorder=1)
    plt.title('Boxplots of $k_{sn}$ by lithology', fontsize=14)
    plt.savefig(DataDirectory+fname_prefix+'_boxplot_lith_ksn.png', dpi=300, transparent=True)
    plt.clf()



    # Do some stats, yo
    # KS test to see if we can distinguish the distributions at a confidence level of p = 0.05
    # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.kstest.html
    # relief
    keys = list(data.keys())
    values = list(data.values())
    k=0
    for i in range(len(keys)-1):
        for j in range(i+1, len(keys)):
            print("KS test between {} and {}".format(keys[i], keys[j]))
            d, p = stats.ks_2samp(values[i], values[j])
            print(d, p)
            k += 1



DataDirectory = '/home/fiona/pCloudDrive/Data_for_papers/river_clusters/Pozo/'
OutDirectory = DataDirectory+'threshold_0/'
fname_prefix = 'Pozo_DTM_basin_208'
raster_name = 'pozo_geol_WGS84_reclass.tif'
label_list = ['Canada shale', 'Upper BV', 'Lower BV', 'SO breccia']
# BoxPlotByCluster(DataDirectory, OutDirectory, fname_prefix,  raster_name, stream_order=1)
#GetLithologyPercentages(DataDirectory,OutDirectory,fname_prefix,raster_name,stream_order=1)
#PlotBasinsWithHillshade(DataDirectory, OutDirectory, fname_prefix, raster_name,stream_order=2)
MakeBoxPlotsKsnLithology(DataDirectory, fname_prefix, raster_name, theta=0.45, label_list=label_list)
