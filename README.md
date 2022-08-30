# TO DO [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6801606.svg)](https://doi.org/10.5281/zenodo.6801606)

# stem_cell_analysis

## Description

These two scripts allow for different analysis in stem cells.

* [differenciation_measurement.py](https://github.com/imcf-shareables/stem_cell_analysis/blob/main/differenciation_measurements.py) measures intensity in both ectoderm and endoderm and count how many cells are positive for either based on a user defined threshold.
* [coexpression_stem_cells.py](https://github.com/imcf-shareables/stem_cell_analysis/blob/main/coexpression_stem_cells.py) measures intensity of both Concanavalin A and ribosomal RNA in a region around the nuclei. This region corresponds to different dilation radius allowing for measurements in both nuclear membrane and cytoplasm.

On top of these files, a [YAML](https://github.com/imcf-shareables/stem_cell_analysis/blob/main/cellpose_env_sr.yaml) file is joined, allowing for reproduction of the Cellpose environment used for the analysis at the [Imaging Core Facility of Basel](https://www.biozentrum.unibas.ch/facilities/technology-platforms/technology-platforms-a-z/overview/unit/imaging-core-facility).

## Requirements

These scripts are made in Python and requires [Fiji](http://fiji.sc/) to be ran. On top of this, multiple update sites need to be activated following [this guide](https://imagej.net/update-sites/#following-an-update-site): 
* 3DImageSuite
* CSBDeep
* IJPB-Plugins
* BIOP
* TrackMate-Cellpose
* TrackMate-StarDist
* CLIJ & CLIJ2

And a [JAR](https://github.com/ome/omero-insight/releases/tag/v5.7.0) from the [OME-team](https://github.com/ome) was used to access files on the IMCF OMERO server.

Once activated, just drag and drop the script in the main Fiji window and click on the RUN button.

As Fiji is operating system independant, this should run on any Windows, Mac and Linux. However, as the scripts use GPU acceleration for processing and Deep Learning algorithms, this will require a GPU with CUDA support to be fully optimal. 

## Run scripts

### Differenciation measurement

#### Input

The script asks for a folder containing a list of files and a specific extension to look for. An upper threshold for the volume of nucleis is used to discard objects that are found by mistakes during segmentation. One threshold for the ectoderm and one for the endoderm are also set in the script to filter cells considered as positive for either. The script also prompts for one or multiple ROIs to be drawn on the image in order to focus the analysis on complete cell clusters.

#### Runtime 

Once initiated, the script uses [StarDist](https://doi.org/10.1007/978-3-030-00934-2_30) ran through [TrackMate](https://www.biorxiv.org/content/10.1101/2021.09.03.458852v2) to find the nuclei in the selected regions. Found nucleis are then filtered based on the volume entered by the user using the [3DImageJSuite](https://academic.oup.com/bioinformatics/article/29/14/1840/231770?login=true) and if the nuclei is touching the X or Y borders, keeping only complete nucleis. Once filtered, the nuclei are dilated by a radius of 2 to which are subtracted the non-dilated ones to only keep a the cytoplasm in close proximity to the nucleis.

The endoderm intensity is then measured in the nuclei region using the corresponding channel, the ectoderm intensity is measured in the cytoplasm region using the corresponding channel. Cells considered positive for the endoderm intensity and cells positive for ectoderm intensity are counted and reported in the output CSV.

#### Output

The script saves two different CSV:
* One reporting the individual measurements of each cells with the volume of the nuclei, the intensities of both ectoderm and endoderm.
* One summing the results with the number of cells and percentage of positive for ectoderm and positive for endoderm.

### Coexpression stem cells

#### Input

The script expects the images to be stored on an OMERO server and prompts the user for both credentials and either a link to all the images or the image IDs separated by comma. A [Cellpose](https://www.nature.com/articles/s41592-020-01018-x) environment is also required for this script and its path is also a required input. A dilation range is expected allowing for measurements to be run in batch for multiple regions of the cell. The script also prompts for one or multiple ROIs to be drawn on the image in order to focus the analysis on complete cell clusters.

#### Runtime 

Once initiated, the script fetches the images one after the other on the OMERO server and will do the analysis in batch. 

The Concanavalin A (ConA) channel is segmented using a median filter with a radius of 2 followed by a threshold using the [Moments algorithm](https://doi.org/10.1111/j.1749-6632.1965.tb11715.x).

The rRNA channel (rRNA) is segmented using a median filter with a radius of 2 followed by a threshold using the [Triangle algorithm](https://doi.org/10.1177/25.7.70454).

The nuclei channel is segmented using [Cellpose](https://www.nature.com/articles/s41592-020-01018-x) ran through [TrackMate](https://www.biorxiv.org/content/10.1101/2021.09.03.458852v2) in the selected regions. Found nucleis are then filtered if touching the X or Y borders of the image using the [3DImageJSuite](https://academic.oup.com/bioinformatics/article/29/14/1840/231770?login=true), keeping only complete nucleis. Once filtered, for each value of the dilation range, nuclei are dilated by that radius using [CLIJ](https://doi.org/10.1038/s41592-019-0650-1) to which we subtract the original nuclei, allowing for measurements to happen in either close proximity to the nuclei or in the whole cytoplasm, depending on the radius. 

The region selected is then combined to the previous segmentations for ConA and rRNA, keeping only detected signal in both channel in the current region of interest. [JACoP](https://doi.org/10.1111/j.1365-2818.2006.01706.x) is then run on these two newly created images and the results stored for output.

#### Output

The script saves two different CSV:
* One per image, containing the image name, the different dilation radiuses, the nuclei intensity as well as the mean intensity in the selected regions for ConA and rRNA
* One for the whole run, containing all the images names, the dilation radiuses, the Pearson and Manders coefficients for each radius.

One TIFF file of the segmented regions per radius can also be saved if the checkbox was ticked.
