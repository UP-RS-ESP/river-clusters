# river-clusters

This repository contains code for clustering rivers in different catchments based on their long profiles. These instructions are for install on a linux machine. For help installing on other platforms, feel free to contact Fiona Clubb (details below).

## Licence

GNU General Public License, Version 3, 29 June 2007

Copyright © 2018 Helmholtz Centre Potsdam GFZ German Research Centre for Geosciences, Potsdam, Germany and University of Kiel, Kiel, Germany.

river-clusters is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. river-clusters is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

## Dependencies

In order to run the clustering code we first of all need to extract all of the river profiles from a catchment using LSDTopoTools.  First of all make sure you have some dependencies installed:
```
sudo apt-get install -y git
sudo apt-get install -y gdal-bin
sudo apt-get install -y python-gdal
sudo apt-get install -y libfftw3-dev
sudo apt-get install -y cmake
```
Then make a directory for LSDTopoTools and clone the LSDTopoTools2 git repository:
```
cd ~
mkdir LSDTopoTools 
cd LSDTopoTools
git clone https://github.com/LSDtopotools/LSDTopoTools2
cd LSDTopoTools2
```
Then install LSDTopoTools2 by running the following within the `LSDTopoTools2` directory:
```
bash lsdtt2_setup.sh
```
This will have installed LSDTopoTools2 onto your machine. 

**NOTE - if you close the terminal you may not have the correct path to LSDTopoTools2 any more. If this is the case, simply run the `lsdtt2_setup.sh` script again in your new terminal before performing the river profile extraction**

## Download and install

The best way to run the code is to use `conda`. Install Anaconda or Miniconda with Python 3.6: https://conda.io/miniconda.html

First of all clone the repository:
```
git clone https://github.com/UP-RS-ESP/river-clusters.git
cd river-clusters
```
Then create a new conda environment using the included `environment.yml` file:
```
conda env create -f environment.yml
```
and activate the new environment by:
```
source activate river-clusters
```
Deactivate the environment at any time with
```
source deactivate river-clusters
```
You can also run the code without conda by ensuring that you have all of the dependencies installed e.g. via pip (see list in the `environment.yml` file.

**NOTE - the conda environment creation can sometimes fail because of a single package. If so, then just remove that package from the list in the `environment.yml` file and then install it afterwards using `conda install <package-name>`**

## Look at the example data

We have provided an example run in the directory `example_data`. This DEM is from a landscape evolution model with a simple contrast in erodibility between the N and S half of the raster. To run the river cluster code using `LSDTopoTools` you need to provide a parameter file. You can find the parameter file for the example data at: `example_data/river_clusters.driver`. It should look like this:

```
# This is a driver file for LSDTopoTools
# Any lines with the # symbol in the first row will be ignored

# File information
dem read extension: bil
dem write extension: bil
read path: /home/clubb/GitHub/river-clusters/example_data/
write path: /home/clubb/GitHub/river-clusters/example_data/
read fname: spatial_K
write fname: spatial_K
threshold_contributing_pixels: 1000

# Parameters for DEM processing
raster_is_filled: false
min_slope_for_fill: 0.001

# What basins do you want?
select_basins_by_order: true
basin_order: 5

# What analyses you want to do
print_all_tributaries_to_csv: true

# Parameters to print data
print_filled_raster: true
print_basin_raster: true
print_junctions_to_csv: true
write hillshade: true
print_channels_to_csv: true
print_stream_order_raster: true
print_junction_index_raster: true
```
**NOTE: You MUST edit the read path and the write path to point to the correct directory on your machine!**

This parameter file should be saved in the *SAME DIRECTORY* as your DEM.

The parameter `basin_order` controls which catchments are selected for clustering.  For the example data, we chose 5th order so that we just get the one basin in the domain. You can adjust this based on your specific requirements.

## Run LSDTopoTools to get the river profiles

Go into the directory with the example data and then run the river profile extraction:
```
cd example_data
lsdtt-river-clusters ./ river_clusters.driver
```
This always follows the format `lsdtt-river-clusters <path-to-data> <name-of-DEM>`

This should have produced the file `DEM_name_all_tribs.csv` along with some other rasters (such as the stream network and a hillshade, for example). So for the example data, your directory should now look like this:
```
$ ls
river_clusters.driver    spatial_K_CN.csv    spatial_K_hs.hdr               spatial_K_SO.bil
spatial_K_all_tribs.csv  spatial_K_Fill.bil  spatial_K_ingestedParam.param  spatial_K_SO.hdr
spatial_K_basins.bil     spatial_K_Fill.hdr  spatial_K_JI.bil
spatial_K_basins.hdr     spatial_K.hdr       spatial_K_JI.hdr
spatial_K.bil            spatial_K_hs.bil    spatial_K_JN.csv

```

## Run the clustering code

To run the clustering code, go back to the main repository and activate the river clustering environment (if using conda):
```
source activate river-clusters
```
Then run:
```
python cluster-river-profiles.py -h
```
which will bring up a help menu with the options.  The standard format is:
```
python cluster-river-profiles.py -dir </path/to/data/folder/> -fname <DEM_name> -so <stream_order_of_choice>
```
For the example data, this command would look like:
```
python cluster-river-profiles.py -dir ./example_data/ -fname spatial_K -so 1
```
## Output

After you have run the python script with the clustering, you should have produced some new data files and plots which you can use to examine the results. Within the main folder `example_data` you should have the following:
* A plot of all the raw profiles: `spatial_K_profiles_upstream.png`
* A plot of all the profiles of the specified stream order: `spatial_K_profiles_SO1.png`
* A longitudinal river profile of the longest channel: `spatial_K_trunk_profile.png`
* A plot of the number of clusters vs. the distance between clusters, which can be used to aid in determining an appropriate number of clusters: `spatial_K_clusters_dist.png`
* A CSV file of the river profiles with the calculated channel gradient for each node: `spatial_K_slopes.csv`
* A CSV file reporting the parameters that you used to run the clustering for reproducibility: `spatial_K_report.csv`
 
The clustering is run at two different threshold levels, which give two different numbers of clusters (see the paper for more details, or contact me).  Therefore there will now be two different directories within `example_data`: `threshold_0` and `threshold_1`. Within each of these directories you should have:
* A csv file with the river profiles and an assigned cluster ID: `spatial_K_clustered_SO1.csv`, where `SO1` means you clustered the first order streams.
* A dendrogram showing the results of the clustering: `spatial_K_dendrogram_SO1.png`
* Two plots showing the profiles coloured by their cluster ID in map view: `spatial_K_elev_clusters_SO1.png` and `spatial_K_hs_clusters_SO1.png`
* The gradient vs. distance profiles separated by cluster. This will be a separate png image for each cluster: e.g. `spatial_K_profiles_SO1_CL1.png` and `spatial_K_profiles_SO1_CL2.png`
* The median gradient vs. distance profiles for all clusters: `spatial_K_profiles_median_SO1.png`
* A boxplot showing the distribution of channel gradient in each cluster: `spatial_K_boxplot_SO1.png`
* Slope-area plots of the channels individually for each cluster with the extracted channel steepness: `spatial_K_SA_median_SO1.png`

## Contact

For more information please contact Fiona Clubb at the University of Potsdam: `clubb@uni-potsdam.de`
