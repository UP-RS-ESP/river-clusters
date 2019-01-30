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
import os
import sys
import pandas as pd
import clustering as cl
import plotting as pl
import raster_plotting as rpl

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
    print("   python cluster-river-profiles.py -h\n")
    print("=======================================================================\n\n ")

if __name__ == '__main__':

    # If there are no arguments, send to the welcome screen
    if (len(sys.argv) < 1):
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()

    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")

    # The options for clustering
    parser.add_argument("-len", "--profile_len", type=int, help="The minimum length of a profile to keep it. Default = 5 nodes.", default=5)
    parser.add_argument("-sw", "--slope_window", type=int, help="The window size for calculating the slope based on a regression through an equal number of nodes upstream and downstream of the node of interest. This is the total number of nodes that are used for calculating the slope. For example, a slope window of 25 would fit a regression through 12 nodes upstream and downstream of the node, plus the node itself. The default is 25 nodes.", default=25)
    parser.add_argument("-m", "--method", type=str, help="The method for clustering, see the scipy linkage docs for more information. The default is 'ward'.", default='ward')
    parser.add_argument("-step", "--step", type=int, help="The regular spacing in metres that you want the profiles to have for the clustering. This should be greater than sqrt(2* DataRes^2).  The default is 2 m which is appropriate for grids with a resolution of 1 m.", default = 2)
    parser.add_argument("-so", "--stream_order", type=int, help="The stream order that you wish to cluster over. Default is 1.", default=1)
    parser.add_argument("-zmax", "--maximum_elevation_for_plotting", type=float, default = 100, help="This is the maximum elevation in the colourbar of the landscape plot.")

    # Options for slope area analysis for comparison
    parser.add_argument("-SA", "--slope_area", type=bool, help='Set to true to make slope-area plots', default=False)

    # Options for plotting catchment metrics. You need to have run an additional LSDTopoTools analysis for this
    parser.add_argument("-CM", "--catchment_metrics", type=bool, default=False, help='Set true to make boxplots of catchment metrics. You need to have run an additional LSDTopoTools analysis for this')

    # Options for raster plotting
    parser.add_argument("-shp", "--shp", type=str, help="Pass a shapefile with the geology for plotting. If nothing is passed then we don't make this plot.", default=None)
    parser.add_argument("-field", "--lith_field", type=str, help="The field name from the shapefile which contains the lithology information", default="geol")
    parser.add_argument("-geol", "--geol_raster", type=str, help="Pass a raster with the geology for plotting.")

    # In case you want to switch the colours. Only works for a two cluster case
    parser.add_argument("-sc", "--switch_colours", type=bool, help="Set to true to switch the colours. Only works for a two cluster case", default=False)

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

    # set min and max of colourbar
    cbar_min_max = [0,args.maximum_elevation_for_plotting]

    # check if the slopes file exists
    slope_file = DataDirectory+args.fname_prefix+'_slopes.csv'
    if os.path.isfile(slope_file):
        df = pd.read_csv(slope_file)
    else:
        # read in the original csv
        df = pd.read_csv(DataDirectory+args.fname_prefix+'_all_tribs.csv')

        # remove profiles with short unique section
        # calculate the slope
        df = cl.CalculateSlope(DataDirectory, args.fname_prefix, df, args.slope_window)
        df.to_csv(DataDirectory+args.fname_prefix+'_slopes.csv', index=False)

    # slope-area plotting if required
    if args.slope_area:
        pl.PlotSlopeArea(DataDirectory, args.fname_prefix)

    pl.PlotTrunkChannel(DataDirectory, args.fname_prefix)

    # get the profiles for the chosen stream order
    new_df = cl.GetProfilesByStreamOrder(DataDirectory, args.fname_prefix, df, args.step, args.slope_window, args.stream_order)
    if args.stream_order > 1:
        new_df = cl.RemoveNonUniqueProfiles(new_df)

    new_df = cl.RemoveProfilesShorterThanThresholdLength(new_df, args.profile_len)
    #
    # do the clustering. We will do this at two threshold levels for the cutoff point.
    thr_levels = [0,1]
    for i in thr_levels:
        print("========================================================")
        print("Running the clustering with threshold level {}".format(i))
        print("========================================================")
        new_dir = DataDirectory+'threshold_{}/'.format(str(i))
        if not os.path.isdir(new_dir):
             os.makedirs(new_dir)
        cl.ClusterProfilesVaryingLength(DataDirectory, new_dir, args.fname_prefix, new_df, args.method, args.stream_order, i)
        if args.switch_colours:
            pl.switch_colours(new_dir, args.fname_prefix, args.stream_order)
        # these functions make some plots for you.
        pl.PlotProfilesByCluster(DataDirectory, new_dir, args.fname_prefix, args.stream_order)
        #rpl.PlotElevationWithClusters(DataDirectory, new_dir, args.fname_prefix, args.stream_order)
        rpl.PlotHillshadewithClusters(DataDirectory, new_dir, args.fname_prefix, args.stream_order)
        if (args.shp == True):
            rpl.PlotLithologyWithClusters(DataDirectory, new_dir, args.fname_prefix, args.stream_order, args.shp, args.lith_field)
        if (args.geol_raster == True):
            rpl.PlotRasterLithologyWithClusters(DataDirectory, new_dir, args.fname_prefix, args.stream_order, args.geol_raster)
        pl.PlotSlopeAreaAllProfiles(DataDirectory, new_dir, args.fname_prefix, args.stream_order, orientation='vertical', nbins=10)
        pl.PlotMedianProfiles(DataDirectory, new_dir, args.fname_prefix, args.stream_order)
        pl.MakeBoxPlotByCluster(DataDirectory, new_dir, args.fname_prefix, args.stream_order)
        if (args.catchment_metrics == True):
            pl.MakeCatchmentMetricsBoxPlot(DataDirectory, new_dir, args.fname_prefix, args.stream_order)

    print('Enjoy your clusters, pal')
