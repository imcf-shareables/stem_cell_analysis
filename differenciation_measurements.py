# ─── SCRIPT PARAMETERS ──────────────────────────────────────────────────────────

#@ File(label="Select the directory with your images", style="directory") src_dir
#@ String(label="Extension for the images to look for", value="czi") filename_filter
#@ Integer(label="Maximum volume for the nuclei", value=3000) max_volume_nuclei
#@ RoiManager rm

# ─── REQUIREMENTS ───────────────────────────────────────────────────────────────

# List of update sites needed for the code
# * 3DImageSuite
# * CSBDeep
# * ImageScience
# * IJPB-Plugins
# * StarDist
# * TrackMate-StarDist

# ─── IMPORTS ────────────────────────────────────────────────────────────────────

import os
import sys
import csv
from itertools import izip

from ij import IJ
from ij.gui import WaitForUserDialog, YesNoCancelDialog
from ij.plugin import Duplicator, RGBStackMerge, ZProjector, ImageCalculator

# Bioformats imports
from loci.plugins import BF, LociExporter
from loci.plugins.in import ImporterOptions
from loci.plugins.out import Exporter

# java imports
from java.lang import Double

# TrackMate imports
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.features import FeatureFilter
from fiji.plugin.trackmate.stardist import StarDistDetectorFactory
from fiji.plugin.trackmate.features.track import TrackDurationAnalyzer
from fiji.plugin.trackmate.action import LabelImgExporter
from fiji.plugin.trackmate.providers import SpotAnalyzerProvider
from fiji.plugin.trackmate.providers import SpotMorphologyAnalyzerProvider

# 3DSuite imports
from mcib3d.geom import Objects3DPopulation
from mcib3d.image3d import ImageInt, ImageHandler

from inra.ijpb.label import LabelImages


# ─── FUNCTIONS ──────────────────────────────────────────────────────────────────
def getFileList(directory, filteringString):
    """Get a list of files with the extension

    Parameters
    ----------
    directory : str
        Path of the files to look at
    filteringString : str
        Extension to look for

    Returns
    -------
    list
        List of files with the extension in the folder
    """

    files = []
    for (dirpath, dirnames, filenames) in os.walk(directory):
    	# if out_dir in dirnames: # Ignore destination directory
            # dirnames.remove(OUT_SUBDIR)
        files.extend(
            os.path.join(dirpath, f) for f in filenames if f.endswith(filteringString)
        )

    return (files)

def BFImport(indivFile):
    """Import using BioFormats

    Parameters
    ----------
    indivFile : str
        Path of the file to open

    Returns
    -------
    imps : ImagePlus
        Image opened via BF
    """
    options = ImporterOptions()
    options.setId(str(indivFile))
    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE)
    options.setOpenAllSeries(True)
    return BF.openImagePlus(options)


def BFExport(imp, savepath):
    """Export using BioFormats

    Parameters
    ----------
    imp : ImagePlus
        ImagePlus of the file to save
    savepath : str
        Path where to save the image

    """
    paramstring = "outfile=[" + savepath + "] windowless=true compression=Uncompressed saveROI=false"


    print('Savepath: ', savepath)
    plugin     = LociExporter()
    plugin.arg = paramstring
    exporter   = Exporter(plugin, imp)
    exporter.run()

def run_tm(implus, channel_number, quality_thresh, intensity_thresh, circularity_thresh, area_thresh, crop_roi=None):
    # sourcery skip: merge-else-if-into-elif, swap-if-else-branches
    """Function to run TrackMate on open data. Has some specific input

    Parameters
    ----------
    implus : ImagePlus
        ImagePlus on which to run image
    channel_number : int
        Number of the channel of interest
    quality_thresh : float
        Value to enter as threshold in the filters
    intensity_thresh : float
        Value to enter as threshold in the filters
    circularity_thresh : float
        Value to enter as threshold in the filters
    area_thresh : float
        Value to enter as threshold in the filters
    crop_roi : ROI, optional
        ROI to crop on the image, by default None
    """

    dims = implus.getDimensions()
    cal = implus.getCalibration()

    if implus.getNSlices() > 1:
        implus.setDimensions(dims[2], dims[4], dims[3])

    # implus2 = implus.duplicate()
    # implus2.setCalibration(implus.getCalibration())

    if crop_roi is not None:
        implus.setRoi(crop_roi)

    model = Model()

    # Send all messages to ImageJ log window.
    model.setLogger(Logger.IJ_LOGGER)

    #------------------------
    # Prepare settings object
    #------------------------

    settings = Settings(implus)

    settings.detectorSettings = {'TARGET_CHANNEL' : channel_number}

    #settings.addAllAnalyzers()
    spotAnalyzerProvider = SpotAnalyzerProvider(1)
    spotMorphologyProvider = SpotMorphologyAnalyzerProvider(1)

    for key in spotAnalyzerProvider.getKeys():
        settings.addSpotAnalyzerFactory( spotAnalyzerProvider.getFactory( key ) )

    for key in spotMorphologyProvider.getKeys():
        settings.addSpotAnalyzerFactory( spotMorphologyProvider.getFactory( key ) )

    #settings.initialSpotFilterValue=quality_thresh

    # Configure detector - We use the Strings for the keys
    settings.detectorFactory = StarDistDetectorFactory()

    # Add the filter on mean intensity
    # Here 'true' takes everything ABOVE the mean_int value
    if (quality_thresh != 0):
        filter1 = FeatureFilter('QUALITY', Double(quality_thresh), True)
        settings.addSpotFilter(filter1)
    if (intensity_thresh != 0):
        filter2 = FeatureFilter('MEAN_INTENSITY_CH' + str(channel_number), Double(intensity_thresh), True)
        settings.addSpotFilter(filter2)
    if (circularity_thresh != 0):
        filter3 = FeatureFilter('CIRCULARITY', Double(circularity_thresh), True)
        settings.addSpotFilter(filter3)
    if (area_thresh != 0):
        filter4 = FeatureFilter('AREA', Double(area_thresh), False)
        settings.addSpotFilter(filter4)


    # Configure tracker
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap()
    settings.addTrackAnalyzer(TrackDurationAnalyzer())
    settings.trackerSettings['LINKING_MAX_DISTANCE'] = 15.0
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 15.0
    settings.trackerSettings['MAX_FRAME_GAP'] = 3

    # print(str(settings))

    #-------------------
    # Instantiate plugin
    #-------------------

    trackmate = TrackMate(model, settings)

    #--------
    # Process
    #--------

    ok = trackmate.checkInput()
    if not ok:
        sys.exit(str(trackmate.getErrorMessage()))

    ok = trackmate.process()
    if not ok:
        if "[SparseLAPTracker] The spot collection is empty." in str(trackmate.getErrorMessage()):
            return IJ.createImage("Untitled", "8-bit black", imp.getWidth(), imp.getHeight(), imp.getNFrames())
        else:
            sys.exit(str(trackmate.getErrorMessage()))

    trackmate.computeSpotFeatures( True )

    # A selection.
    sm = SelectionModel( model )

    # Read the default display settings.
    # ds = DisplaySettingsIO.readUserDefault()
    # displayer =  HyperStackDisplayer(model, sm, implus, ds)
    # displayer.render()
    # displayer.refresh()

    exportSpotsAsDots = False
    exportTracksOnly = False
    # implus2.close()
    temp= LabelImgExporter.createLabelImagePlus( trackmate, exportSpotsAsDots, exportTracksOnly )
    temp.setDimensions(dims[2], dims[3], dims[4])
    implus.setDimensions(dims[2], dims[3], dims[4])
    return temp

# ─── VARIABLES ──────────────────────────────────────────────────────────────────

filter_objects_touching_z = False

# Channel order
endoderm_chnl = 1
ectoderm_chnl = 2
dapi_chnl     = 3

endo_thresh = 4000
ecto_thresh = 2000

# Variable for trackmate, is ran on 2D labels
quality_thr     = 0
intensity_thr   = 0
circularity_thr = 0
area_thr        = 0

# ─── MAIN CODE ──────────────────────────────────────────────────────────────────
IJ.log("\\Clear")
IJ.log("Script starting")

# Retrieve list of files
src_dir = str(src_dir)
files = getFileList(src_dir, filename_filter)
rm.reset()

# # If the list of files is not empty
if files:

    # For each file finishing with the filtered string
    for file in sorted(files):

        # Get info for the files
        folder   = os.path.dirname(file)
        basename = os.path.basename(file)
        basename = os.path.splitext(basename)[0]



        # Import the file with BioFormats
        IJ.log("Currently opening " + basename + "...")

        imps = BFImport(str(file))
        for imp in imps:

            rm.reset()
            basename = imp.getTitle()

            out_roi = os.path.join(folder, basename + "_ROI.zip")
            out_single_roi = os.path.join(folder, basename + "_ROI.roi")
            out_3D_rois = os.path.join(folder, basename + "_3DROIs.zip")
            out_label = os.path.join(folder, basename + "_nuclei_labels.tif")
            out_csv = os.path.join(folder, basename + "_results.csv")
            out_csv2 = os.path.join(folder, basename + "_image_results.csv")

            if (os.path.exists(out_roi) or os.path.exists(out_single_roi)):
            # query_roi = YesNoCancelDialog(None, "ROI file exists","The file" +
            #                         " containing a drawn ROI exists," +
            #                         " this means that the script" +
            #                         " has already been ran. Do you want to" +
            #                         " run it again on that file ?",
            #                         "Reuse", "Draw new one")
            # if query_roi.yesPressed():
                if os.path.exists(out_single_roi):
                    rm.runCommand("Open", out_single_roi)
                elif os.path.exists(out_roi):
                    rm.runCommand("Open", out_roi)
                # roi = rm.getRoi(0)
                    # imp.setRoi(roi)
                # if query_roi.cancelPressed():
                #     sys.exit("Cancel")

            cal    = imp.getCalibration()
            dims   = imp.getDimensions()
            unit   = cal.getUnits()

            count = 0
            while(rm.getCount() == 0):
                if (count == 0):
                    imp.show()

                WaitForUserDialog(
                    "Draw the region(s) of interest, press T to add to ROI Manager and press OK").show()
                if count == 5:
                    imp.close()
                    sys.exit("Too many clicks without ROI")
                else:
                    count += 1

            if (not os.path.exists(out_roi) or not os.path.exists(out_single_roi)):
                rm.runCommand("Save", out_roi)

            imp.hide()
            IJ.run(imp, "Select None", "")

            endoderm_imp = Duplicator().run(imp,
                                            endoderm_chnl, endoderm_chnl,
                                            1, imp.getNSlices(),
                                            1, imp.getNFrames())
            ectoderm_imp = Duplicator().run(imp,
                                            ectoderm_chnl, ectoderm_chnl,
                                            1, imp.getNSlices(),
                                            1, imp.getNFrames())

            for roi_nb in range(rm.getCount()):
                current_roi = rm.getRoi(roi_nb)
                IJ.log("Current working on ROI " + current_roi.getName())

                IJ.run(imp, "Select None", "")

                current_label_imp = run_tm(imp, dapi_chnl, quality_thr, intensity_thr,
                                        circularity_thr, area_thr, current_roi)
                current_label_imp.setTitle("Current_label_imp")

                if roi_nb == 0:
                    label_imp = Duplicator().run(current_label_imp,
                                                    1, 1,
                                                    1, current_label_imp.getNSlices(),
                                                    1, current_label_imp.getNFrames())
                    label_imp.setTitle("Label_imp")
                    current_label_imp.close()
                else:

                    mip_temp = ZProjector().run(label_imp,"max")
                    current_objects_number = mip_temp.getStatistics().max
                    current_label_imp2 = Duplicator().run(current_label_imp,
                                                    1, 1,
                                                    1, current_label_imp.getNSlices(),
                                                    1, current_label_imp.getNFrames())
                    IJ.run(current_label_imp2, "Add...",
                            "value=" + str(current_objects_number) + " stack")
                    IJ.run(current_label_imp2, "Replace/Remove Label(s)",
                            "label(s)=" + str(current_objects_number) +" final=0")

                    ImageCalculator.run(label_imp, current_label_imp2,
                                                    "Add stack")

                    current_label_imp.close()
                    current_label_imp2.close()
                    mip_temp.close()

            img_nuclei              = ImageInt.wrap(label_imp)
            pop_nuclei              = Objects3DPopulation(img_nuclei)
            nb_nuclei               = pop_nuclei.getNbObjects()
            nucleis_obj_to_remove   = []
            nucleis_index_to_remove = []

            # Loop through objects to remove too big ones
            for i in range(0, nb_nuclei):
                obj = pop_nuclei.getObject(i)

                if obj.getVolumeUnit() > max_volume_nuclei:
                    nucleis_obj_to_remove.append(obj)
                    nucleis_index_to_remove.append(obj.getValue())
                elif obj.touchBorders(img_nuclei, filter_objects_touching_z):
                    nucleis_obj_to_remove.append(obj)
                    nucleis_index_to_remove.append(obj.getValue())

            for obj in nucleis_obj_to_remove:
                pop_nuclei.removeObject(obj)

            for i in nucleis_index_to_remove:
                IJ.run(label_imp, "Replace/Remove Label(s)",
                            "label(s)=" + str(i) +" final=0")

            dilated_label_imp = Duplicator().run(label_imp,
                                1, 1,
                                1, label_imp.getNSlices(),
                                1, label_imp.getNFrames())

            # dilated_label_imp.show()
            # IJ.run(dilated_label_imp, "Dilate labels", "radius=2")
            dilated_label_imp = LabelImages.dilateLabels(dilated_label_imp, 2)
            ImageCalculator.run(dilated_label_imp, label_imp,
                            "Subtract stack")

            img_cyto = ImageInt.wrap(dilated_label_imp)
            pop_cyto = Objects3DPopulation(img_cyto)
            nb_cyto  = pop_cyto.getNbObjects()


            IH_endoderm = ImageHandler.wrap(endoderm_imp)
            IH_ectoderm = ImageHandler.wrap(ectoderm_imp)

            obj_index = []
            obj_volume = []
            endoderm_intensity = []
            ectoderm_intensity = []

            pos_endo = 0
            pos_ecto = 0

            for i in range(0, nb_cyto):
                nuclei_obj = pop_nuclei.getObject(i)
                cyto_obj = pop_cyto.getObject(i)

                obj_index.append(nuclei_obj.getValue())
                obj_volume.append(nuclei_obj.getVolumeUnit())
                endoderm_intensity.append(nuclei_obj.getPixMeanValue(IH_endoderm))
                ectoderm_intensity.append(cyto_obj.getPixMeanValue(IH_ectoderm))

                if cyto_obj.getPixMeanValue(IH_ectoderm) > ecto_thresh:
                    pos_ecto+=1
                if nuclei_obj.getPixMeanValue(IH_endoderm) > endo_thresh:
                    pos_endo+=1


            # IJ.log("Saving 3D ROIs...")
            # pop_nuclei.saveObjects(out_3D_rois)
            IJ.log("Saving label image...")
            BFExport(label_imp, out_label)

            with open(out_csv, 'wb') as f:
                writer = csv.writer(f)
                writer.writerow(
                    ["Object ID, Volume (" + unit + " cube), Mean intensity endoderm, " +
                    "Mean intensity ectoderm"])
                writer.writerows(izip(obj_index, obj_volume, endoderm_intensity,
                    ectoderm_intensity))
            f.close()

            with open(out_csv2, 'wb') as f:
                writer = csv.writer(f)
                writer.writerow(
                    ["Total number of objects", "% Positive for endo", "% Positive for ecto"]
                )
                writer.writerow([nb_cyto,
                                 round((float(pos_endo) / float(nb_cyto)) * 100,2),
                                 round((float(pos_ecto) / float(nb_cyto)) * 100,2)])
            f.close()

        imp.close()
        label_imp.close()
        dilated_label_imp.close()
        endoderm_imp.close()
        ectoderm_imp.close()



    IJ.log("Script DONE ! ")