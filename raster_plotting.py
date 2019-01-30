#---------------------------------------------------------------------#
# Raster plotting functions for clustering analysis
# Developed by Fiona Clubb
#              Bodo Bookhagen
#              Aljoscha Rheinwalt
# University of Potsdam
#---------------------------------------------------------------------#

from LSDPlottingTools import LSDMap_VectorTools as VT
from LSDPlottingTools import LSDMap_GDALIO as IO
from LSDPlottingTools import LSDMap_BasicManipulation as BM
from LSDMapFigure import PlottingRaster
from LSDMapFigure.PlottingRaster import MapFigure
import pandas as pd
from shapely.geometry import shape, Polygon
from descartes.patch import PolygonPatch
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.cm as cm
from matplotlib import rcParams
import LSDPlottingTools as LSDP

# Set up fonts for plots
label_size = 12
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = label_size


def PlotElevationWithClusters(DataDirectory, OutDirectory, fname_prefix, stream_order=1, cbar_loc='right', custom_cbar_min_max = []):
    """
    Make a plot of the raster with the channels coloured by the cluster
    value. Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

    Args:
        stream_order: the stream order of the profiles that you are analysing
        cbar_loc: location of the colourbar, can be right, top, bottom, left, or none.
        custom_cbar_min_max: list of [min, max] to recast the raster to for display.

    Author: FJC
    """
    print("I'm plotting the elevation with channels coloured by cluster")
    df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
    cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))


    # set figure sizes based on format
    fig_width_inches = 4.92126

    # some raster names
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HSName = fname_prefix+'_hs'+raster_ext

    if not os.path.isfile(DataDirectory+HSName):
        # make a hillshade
        BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

    MF = MapFigure(HSName, DataDirectory,coord_type="UTM",colourbar_location = cbar_loc)
    MF.add_drape_image(BackgroundRasterName,DataDirectory,colourmap = 'gray', alpha=0.5, colorbarlabel = "Elevation (m)",colour_min_max = custom_cbar_min_max)

    # # create the map figure
    # MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM",colour_min_max = custom_cbar_min_max)

    clusters = cluster_df.cluster_id.unique()
    for cl in clusters:
        # plot the whole channel network in black
        ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
        MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='w',manual_size=1, zorder=2, alpha=1)
        # plot the clustered profiles in the correct colour
        this_df = cluster_df[cluster_df.cluster_id == cl]
        this_colour = str(this_df.colour.unique()[0])
        ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
        MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=2)

    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = OutDirectory+fname_prefix+'_elev_clusters_SO{}.png'.format(stream_order), FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False, axis_style='Thin', transparent=True) # Save the figure

def PlotHillshadewithClusters(DataDirectory, OutDirectory, fname_prefix,stream_order=1):
        """
        Make a hillshade of the raster with the channels coloured by the cluster
        value. Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

        Args:
            stream_order: the stream order of the profiles that you are analysing

        Author: FJC
        """
        import LSDPlottingTools as LSDP
        from LSDMapFigure.PlottingRaster import MapFigure

        print("I'm plotting a shaded relief map with the channels coloured by cluster")

        df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
        cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))


        # set figure sizes based on format
        fig_width_inches = 8

        # some raster names
        raster_ext = '.bil'
        BackgroundRasterName = fname_prefix+raster_ext
        HSName = fname_prefix+'_hs'+raster_ext

        if not os.path.isfile(DataDirectory+HSName):
            # make a hillshade
            BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

        # create the map figure
        MF = MapFigure(HSName, DataDirectory,coord_type="UTM")

        clusters = cluster_df.cluster_id.unique()
        for cl in clusters:
            # plot the whole channel network in black
            ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
            MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='0.9',manual_size=2, zorder=1, alpha=0.5)
            # plot the clustered profiles in the correct colour
            this_df = cluster_df[cluster_df.cluster_id == cl]
            this_colour = str(this_df.colour.unique()[0])
            ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
            MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=3)

        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = OutDirectory+fname_prefix+'_hs_clusters_SO{}.png'.format(stream_order), FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False, transparent=True) # Save the figure

def PlotLithologyWithClusters(DataDirectory, OutDirectory, fname_prefix, stream_order=1, shapefile_name = 'geol.shp', geol_field = 'geol'):
        """
        Make a hillshade of the raster with the channels coloured by the cluster
        value. Rasterise a geology shapefile and drape on the top.
        Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

        Args:
            stream_order: the stream order of the profiles that you are analysing
            shapefile_name: name of the lithology shapefile
            geol_field: the field of the shapefile that has the lithology information

        Author: FJC
        """
        import LSDPlottingTools as LSDP
        from LSDMapFigure.PlottingRaster import MapFigure

        df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
        cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))


        # set figure sizes based on format
        fig_width_inches = 4.92126

        # some raster names
        raster_ext = '.bil'
        BackgroundRasterName = fname_prefix+raster_ext
        HSName = fname_prefix+'_hs'+raster_ext

        if not os.path.isfile(DataDirectory+HSName):
            print("Making a hillshade for you")
            # make a hillshade
            BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

        # create the map figure
        MF = MapFigure(HSName, DataDirectory,coord_type="UTM")
        res = IO.GetUTMMaxMin(DataDirectory+BackgroundRasterName)[0]

        #rasterise the shapefile
        new_shp, geol_dict = VT.geologic_maps_modify_shapefile(DataDirectory+ shapefile_name, geol_field)
        lith_raster = VT.Rasterize_geologic_maps_pythonic(new_shp, res, geol_field)
        MF.add_drape_image(lith_raster, "", colourmap=plt.cm.jet, alpha=0.4, show_colourbar = False, discrete_cmap=True, cbar_type=int,mask_value=0)

        clusters = cluster_df.cluster_id.unique()
        for cl in clusters:
            # plot the whole channel network in black
            ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
            MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='0.9',manual_size=2, zorder=1, alpha=0.5)
            # plot the clustered profiles in the correct colour
            this_df = cluster_df[cluster_df.cluster_id == cl]
            this_colour = str(this_df.colour.unique()[0])
            ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
            MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=3)

        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = OutDirectory+fname_prefix+'_lith_clusters_SO{}.png'.format(stream_order), FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False, transparent=True) # Save the figure


def PlotRasterLithologyWithClusters(DataDirectory, OutDirectory, fname_prefix, stream_order=1, geol_raster = 'geol'):
        """
        Make a hillshade of the raster with the channels coloured by the cluster
        value. Rasterise a geology shapefile and drape on the top.
        Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

        Args:
            stream_order: the stream order of the profiles that you are analysing
            shapefile_name: name of the lithology shapefile
            geol_field: the field of the shapefile that has the lithology information

        Author: FJC
        """
        import LSDPlottingTools as LSDP
        from LSDMapFigure.PlottingRaster import MapFigure

        df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
        cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))


        # set figure sizes based on format
        fig_width_inches = 8

        # some raster names
        raster_ext = '.bil'
        BackgroundRasterName = fname_prefix+raster_ext
        HSName = fname_prefix+'_hs'+raster_ext

        if not os.path.isfile(DataDirectory+HSName):
            # make a hillshade
            BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

        # create the map figure
        MF = MapFigure(HSName, DataDirectory,coord_type="UTM")

        #geology
        LithName = geol_raster
        print("The geology raster is"+LithName)
        MF.add_drape_image(LithName, DataDirectory, colourmap=plt.cm.jet, alpha=0.5, show_colourbar = False, discrete_cmap=True, cbar_type=int,mask_value=0)

        clusters = cluster_df.cluster_id.unique()
        for cl in clusters:
            # plot the whole channel network in black
            ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
            MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='white',manual_size=1.5, zorder=2, alpha=0.5)
            # plot the clustered profiles in the correct colour
            this_df = cluster_df[cluster_df.cluster_id == cl]
            this_colour = str(this_df.colour.unique()[0])
            ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
            MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=2.5)

        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = OutDirectory+fname_prefix+'_lith_clusters_SO{}.png'.format(stream_order), FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False, transparent=True) # Save the figure

def PlotRasterLithology(DataDirectory, fname_prefix, geol_raster = 'geol'):
        """
        Make a hillshade and rasterise a geology shapefile and drape on the top.
        Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

        Args:
            stream_order: the stream order of the profiles that you are analysing
            shapefile_name: name of the lithology shapefile
            geol_field: the field of the shapefile that has the lithology information

        Author: FJC
        """
        import LSDPlottingTools as LSDP
        from LSDMapFigure.PlottingRaster import MapFigure

        # set figure sizes based on format
        fig_width_inches = 8

        # some raster names
        raster_ext = '.bil'
        BackgroundRasterName = fname_prefix+raster_ext
        HSName = fname_prefix+'_hs'+raster_ext

        if not os.path.isfile(DataDirectory+HSName):
            # make a hillshade
            BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

        # create the map figure
        MF = MapFigure(HSName, DataDirectory,coord_type="UTM")

        #geology
        LithName = geol_raster
        print("The geology raster is"+LithName)
        MF.add_drape_image(LithName, DataDirectory, colourmap=plt.cm.jet, alpha=0.5, show_colourbar = False, discrete_cmap=True, cbar_type=int,mask_value=0)

        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = DataDirectory+fname_prefix+'_lith.png', FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False, transparent=True) # Save the figure

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
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,5:100])

    # plot the raster
    hs_raster = IO.ReadRasterArrayBlocks(DataDirectory+fname_prefix+'_hs.bil')
    extent = IO.GetRasterExtent(DataDirectory+fname_prefix+'_hs.bil')
    # hs_raster = IO.ReadRasterArrayBlocks(DataDirectory+'Pozo_DTM_basin_208_hs.bil')
    # extent = IO.GetRasterExtent(DataDirectory+'Pozo_DTM_basin_208_hs.bil')
    plt.imshow(hs_raster, cmap=cm.gray, extent=extent)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel('Easting (m)')
    plt.ylabel('Northing (m)')

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

    plt.savefig(OutDirectory+fname_prefix+'_hs_basins_SO{}.png'.format(stream_order), FigFormat='png', dpi=500, transparent=True)

def PlotKsnFromSlopeArea(DataDirectory, fname_prefix, theta=0.45, cbar_loc='right', custom_cbar_min_max = []):
    """
    Make a plot of the slope area data with a fixed concavity
    """
    import numpy as np

    print("Calculating ksn...")
    # read the csv and get some info
    df = pd.read_csv(DataDirectory+fname_prefix+"_slopes.csv")

    # now force a fit of ks based on this concavity
    area = df['drainage_area'].values
    slope = df['slope'].values

    ksn = slope/(area**(-theta))
    df['ksn'] = ksn
    print(np.max(ksn), np.min(ksn))
    print(ksn)

    df.to_csv(DataDirectory+fname_prefix+'_ksn.csv', index=False)

    # set figure sizes based on format
    fig_width_inches = 6

    # some raster names
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HSName = fname_prefix+'_hs'+raster_ext

    if not os.path.isfile(DataDirectory+HSName):
        # make a hillshade
        BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

    # create the map figure
    MF = MapFigure(HSName, DataDirectory,coord_type="UTM", cbar_loc='right', font_size=16)
    #MF.add_drape_image(BackgroundRasterName,DataDirectory,colourmap = 'gray', alpha=0.5)

    # add the ksn data
    ChannelPoints = LSDP.LSDMap_PointData(df, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelPoints, this_colourmap='viridis', column_for_plotting='ksn', colour_log=True,zorder=100, show_colourbar = True, colourbar_location='bottom', font_size=24, colorbarlabel = "log$_{10}(k_{sn})", manual_size = 1.5, colour_manual_scale = [0,1.8])
    #plt.show()

    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = DataDirectory+fname_prefix+'_ksn.png', FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False, axis_style='Thin', transparent=True) # Save the figure
    plt.clf()


if __name__ == '__main__':

    DataDirectory = '/home/clubb/OneDrive/river_clusters/Pozo/'
    fname_prefix  = 'Pozo_DTM'
    stream_order = 1
    shp = 'pozo_geol_WGS84.shp'
    field = 'Lithology'
    #PlotLithologyWithClusters(DataDirectory, fname_prefix, stream_order, shp, field)
    #PlotRasterLithologyWithClusters(DataDirectory, fname_prefix, stream_order, geol_raster='pozo_geol_WGS84_new_padded')
    PlotBasinsWithHillshade(DataDirectory, DataDirectory+'threshold_0/', fname_prefix, stream_order)
