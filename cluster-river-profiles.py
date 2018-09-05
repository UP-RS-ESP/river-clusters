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
    print("   python plot_river_clusters.py -h\n")
    print("=======================================================================\n\n ")

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
    parser.add_argument("-len", "--profile_len", type=int, help="The minimum length of a profile to keep it. Default = 5 nodes.", default=5)
    parser.add_argument("-sw", "--slope_window", type=int, help="The window size for calculating the slope based on a regression through an equal number of nodes upstream and downstream of the node of interest. This is the total number of nodes that are used for calculating the slope. For example, a slope window of 25 would fit a regression through 12 nodes upstream and downstream of the node, plus the node itself. The default is 25 nodes.", default=25)
    parser.add_argument("-m", "--method", type=str, help="The method for clustering, see the scipy linkage docs for more information. The default is 'ward'.", default='ward')
    parser.add_argument("-c", "--min_corr", type=float, help="The minimum correlation for defining the clusters. Use a smaller number to get less clusters, and a bigger number to get more clusters (from 0 = no correlation, to 1 = perfect correlation). The default is 0.5. DEPRECATED - now we calculate the threshold statistically.", default=0.5)
    parser.add_argument("-step", "--step", type=int, help="The regular spacing in metres that you want the profiles to have for the clustering. This should be greater than sqrt(2* DataRes^2).  The default is 2 m which is appropriate for grids with a resolution of 1 m.", default = 2)
    parser.add_argument("-so", "--stream_order", type=int, help="The stream order that you wish to cluster over. Default is 1.", default=1)

    # Options for raster plotting
    parser.add_argument("-shp", "--shp", type=str, help="Pass a shapefile with the geology for plotting. If nothing is passed then we don't make this plot.", default=None)
    parser.add_argument("-field", "--lith_field", type=str, help="The field name from the shapefile which contains the lithology information", default="geol")

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

        # remove profiles with short unique section
        # calculate the slope
        df = cl.CalculateSlope(DataDirectory, args.fname_prefix, df, args.slope_window)
        df.to_csv(DataDirectory+args.fname_prefix+'_slopes.csv', index=False)

    # get the profiles for the chosen stream order
    new_df = cl.GetProfilesByStreamOrder(DataDirectory, args.fname_prefix, df, args.step, args.slope_window, args.stream_order)
    if args.stream_order > 1:
        new_df = cl.RemoveNonUniqueProfiles(new_df)

    new_df = cl.RemoveProfilesShorterThanThresholdLength(new_df, args.profile_len)

    # do the clustering
    cl.ClusterProfilesVaryingLength(DataDirectory, args.fname_prefix, new_df, args.method, args.stream_order)
    pl.PlotProfilesByCluster(DataDirectory, args.fname_prefix, args.stream_order)
    # # #
    # # #PlotMedianProfiles()
    rpl.PlotElevationWithClusters(DataDirectory, args.fname_prefix, args.stream_order)
    if args.shp:
        rpl.PlotLithologyWithClusters(DataDirectory, args.fname_prefix, args.stream_order, args.shp, args.lith_field)
    pl.PlotSlopeArea(DataDirectory, args.fname_prefix, args.stream_order)
    pl.PlotTrunkChannel(DataDirectory, args.fname_prefix)
    #PlotLongitudinalProfiles()
    #MakeShadedSlopeMap()

    print('Enjoy your clusters, pal')

    #PlotUniqueStreamsWithLength()