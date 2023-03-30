# ─── SCRIPT PARAMETERS ──────────────────────────────────────────────────────────

#@ String(label="Username", description="please enter your username") USERNAME
#@ String(label="Password", description="please enter your password", style="password") PASSWORD
#@ String(label="Info about file(s)", description="Link got from OMERO, or image IDs separated by commas") OMERO_link
#@ File(label="Path for storage", style="directory", description="Where to store the results") destination
#@ File(label="Path to cellpose environment", style="directory", description="Path to conda env") cellpose_env
#@ Integer(label="Min dilation radius for the nuclei", value=2) min_dilation_radius
#@ Integer(label="Max dilation radius for the nuclei", value=5) max_dilation_radius
#@ Boolean(label="Save images ?", value=True) save_tiffs
#@ RoiManager rm

# ─── REQUIREMENTS ───────────────────────────────────────────────────────────────

# List of update sites needed for the code
# * CSBDeep
# * IJPB-Plugins
# * BIOP
# * TrackMate-Cellpose
# * OMERO-ij jar from website
# * CLIJ & CLIJ2

# ─── IMPORTS ────────────────────────────────────────────────────────────────────

import os
import sys
import csv
from itertools import izip

from ij import IJ, Prefs
from ij.gui import WaitForUserDialog
from ij.plugin import Duplicator, ZProjector, ImageCalculator

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
from fiji.plugin.trackmate.tracking.jaqaman import LAPUtils, SparseLAPTrackerFactory
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel
from fiji.plugin.trackmate.action import LabelImgExporter


# 3DSuite imports
from mcib3d.geom import Objects3DPopulation
from mcib3d.image3d import ImageInt, ImageHandler

# Omero Dependencies
from omero.gateway import Gateway
from omero.gateway import LoginCredentials
from omero.log import SimpleLogger

from net.haesleinhuepf.clij2 import CLIJ2

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

def run_tm(implus, channel_number, quality_thresh, intensity_thresh, circularity_thresh, area_thresh, cellpose_env, crop_roi=None):
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
    cellpose_env : str
        Path to the cellpose environment
    crop_roi : ROI, optional
        ROI to crop on the image, by default None
    """

    dims = implus.getDimensions()
    cal = implus.getCalibration()

    if implus.getNSlices() > 1:
        implus.setDimensions(dims[2], dims[4], dims[3])

    # implus2 = Duplicator().run(implus, nuclei_chnl, nuclei_chnl, 1, 1, 1, implus.getNFrames())
    # implus2.setCalibration(implus.getCalibration())

    if crop_roi is not None:
        implus.setRoi(crop_roi)

    model = Model()

    # # Send all messages to ImageJ log window.
    model.setLogger(Logger.IJTOOLBAR_LOGGER)

    #------------------------
    # Prepare settings object
    #------------------------

    settings = Settings(implus)


    # Configure detector - We use the Strings for the keys
    settings.detectorFactory = CellposeDetectorFactory()

    settings.detectorSettings['TARGET_CHANNEL'] = nuclei_chnl
    settings.detectorSettings['OPTIONAL_CHANNEL_2'] = 0
    settings.detectorSettings['CELLPOSE_PYTHON_FILEPATH'] = os.path.join(cellpose_env, 'python.exe')
    settings.detectorSettings['CELLPOSE_MODEL_FILEPATH'] = os.path.join(os.environ['USERPROFILE'], '.cellpose', 'models')
    settings.detectorSettings['CELLPOSE_MODEL'] = PretrainedModel.CYTO2
    settings.detectorSettings['CELL_DIAMETER'] = 11.0
    settings.detectorSettings['USE_GPU'] = True
    settings.detectorSettings['SIMPLIFY_CONTOURS'] = True

    #settings.addAllAnalyzers()
    # spotAnalyzerProvider = SpotAnalyzerProvider(1)
    # spotMorphologyProvider = SpotMorphologyAnalyzerProvider(1)

    # for key in spotAnalyzerProvider.getKeys():
    #     settings.addSpotAnalyzerFactory( spotAnalyzerProvider.getFactory( key ) )

    # for key in spotMorphologyProvider.getKeys():
    #     settings.addSpotAnalyzerFactory( spotMorphologyProvider.getFactory( key ) )

    # #settings.initialSpotFilterValue=quality_thresh

    # # Add the filter on mean intensity
    # # Here 'true' takes everything ABOVE the mean_int value
    # if (quality_thresh != 0):
    #     filter1 = FeatureFilter('QUALITY', Double(quality_thresh), True)
    #     settings.addSpotFilter(filter1)
    # if (intensity_thresh != 0):
    #     filter2 = FeatureFilter('MEAN_INTENSITY_CH' + str(channel_number), Double(intensity_thresh), True)
    #     settings.addSpotFilter(filter2)
    # if (circularity_thresh != 0):
    #     filter3 = FeatureFilter('CIRCULARITY', Double(circularity_thresh), True)
    #     settings.addSpotFilter(filter3)
    # if (area_thresh != 0):
    #     filter4 = FeatureFilter('AREA', Double(area_thresh), False)
    #     settings.addSpotFilter(filter4)


    # Configure tracker
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
    # settings.addTrackAnalyzer(TrackDurationAnalyzer())
    settings.trackerSettings['LINKING_MAX_DISTANCE'] = 15.0
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 15.0
    settings.trackerSettings['MAX_FRAME_GAP'] = 3
    settings.initialSpotFilterValue = -1.

    # print(str(settings))

    #-------------------
    # Instantiate plugin
    #-------------------

    trackmate = TrackMate(model, settings)
    trackmate.computeSpotFeatures( True )
    trackmate.computeTrackFeatures( True )

    #--------
    # Process
    #--------

    ok = trackmate.checkInput()
    if not ok:
        sys.exit(str(trackmate.getErrorMessage()))
        return

    ok = trackmate.process()
    if not ok:
        if "[SparseLAPTracker] The spot collection is empty." in str(trackmate.getErrorMessage()):
            return IJ.createImage("Untitled", "8-bit black", imp.getWidth(), imp.getHeight(), imp.getNFrames())
        else:
            sys.exit(str(trackmate.getErrorMessage()))
            return


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
    temp= LabelImgExporter.createLabelImagePlus( trackmate, exportSpotsAsDots, exportTracksOnly, False )
    temp.setDimensions(dims[2], dims[3], dims[4])
    implus.setDimensions(dims[2], dims[3], dims[4])
    return temp

def omero_connect():
    """Connect to OMERO using the credentials entered

    Returns
    -------
    gateway
        OMERO gateway
    """
    # Omero Connect with credentials and simpleLogger
    cred = LoginCredentials()
    cred.getServer().setHostname(HOST)
    cred.getServer().setPort(PORT)
    cred.getUser().setUsername(USERNAME.strip())
    cred.getUser().setPassword(PASSWORD.strip())
    simpleLogger = SimpleLogger()
    gateway = Gateway(simpleLogger)
    gateway.connect(cred)
    return gateway

def parse_url(omero_str):
    """Parse an OMERO URL with one or multiple images selected

    Parameters
    ----------
    omero_str : str
        String which is either the link gotten from OMERO or image IDs separated by commas

    Returns
    -------
    str[]
        List of all the images IDs parsed from the string
    """
    if omero_str.startswith("https"):
        image_ids = omero_str.split("image-")
        image_ids.pop(0)
        image_ids = [s.split('%')[0].replace("|", "") for s in image_ids]
    else:
        image_ids = omero_str.split(",")
    return image_ids

def progress_bar(progress, total, line_number, prefix=''):
    """Progress bar for the IJ log window

    Parameters
    ----------
    progress : int
        Current step of the loop
    total : int
        Total number of steps for the loop
    line_number : int
        Number of the line to be updated
    prefix : str, optional
        Text to use before the progress bar, by default ''
    """

    size = 30
    x    = int(size*progress/total)
    IJ.log("\\Update%i:%s\t[%s%s] %i/%i\r" % (line_number, prefix, "#"*x, "."*(size-x), progress, total))

def openImagePlus(HOST, USERNAME, PASSWORD, groupId, imageId):
    """Open an ImagePlus from an OMERO server

    Parameters
    ----------
    HOST : str
        Adress of your OMERO server
    USERNAME : str
        Username to use in OMERO
    PASSWORD : str
        Password
    groupId : double
        OMERO group ID
    imageId : int
        ID of the image to open
    """
    stackview = "viewhyperstack=true stackorder=XYCZT "
    datasetorg = "groupfiles=false swapdimensions=false openallseries=false concatenate=false stitchtiles=false"
    coloropt = "colormode=Default autoscale=true"
    metadataview = "showmetadata=false showomexml=false showrois=true setroismode=roimanager"
    memorymanage = "virtual=false specifyranges=false setcrop=false"
    split = " splitchannels=false splitfocalplanes=false splittimepoints=false"
    other = "windowless=true"
    options = ("location=[OMERO] open=[omero:server=%s\nuser=%s\npass=%s\ngroupID=%s\niid=%s] %s %s %s %s %s %s %s " %
               (HOST, USERNAME, PASSWORD, groupId, imageId, stackview, datasetorg, coloropt, metadataview, memorymanage, split, other))
    IJ.runPlugIn("loci.plugins.LociImporter", options)

def dilate_labels_on_gpu(clij2_instance, label_image, dilation_radius):

    #clij2 = CLIJ2.getInstance()
    src = clij2_instance.push(label_image)
    dst = clij2_instance.create(src)

    clij2_instance.dilateLabels(src, dst, dilation_radius)
    return clij2_instance.pull(dst)


# ─── VARIABLES ──────────────────────────────────────────────────────────────────

# OMERO server info
HOST    = "omero.server.address"
PORT    = 4064
groupId = "-1"

filter_objects_touching_z = False

# Channel order
nuclei_chnl = 3
conA_chnl   = 1
rrna_chnl   = 2


# Variable for trackmate, is ran on 2D labels
quality_thr     = 0
intensity_thr   = 0
circularity_thr = 0
area_thr        = 0

# ─── MAIN CODE ──────────────────────────────────────────────────────────────────
IJ.log("\\Clear")
IJ.log("Script starting")

clij2_instance = CLIJ2.getInstance()

# Retrieve list of files
destination = destination.toString()
cellpose_env = cellpose_env.toString()
# print(type(cellpose_env))
# files = getFileList(destination, filename_filter)
rm.reset()

try:

    gateway = omero_connect()

    image_ids_array = parse_url(OMERO_link)
    image_ids_array.sort()

    image_title_array = []
    dilation_array    = []
    pearson_array     = []
    m1_array          = []
    m2_array          = []

    coexpr_csv = os.path.join(destination,
                              "Results_coexpr.csv")

    for image_index, image_id in enumerate(image_ids_array):
        IJ.log("\\Clear")
        clij2_instance.clear()

        progress_bar(image_index + 1, len(image_ids_array), 2, "Processing : ")


# # # If the list of files is not empty
# if files:
#     # For each file finishing with the filtered string
#     for file in sorted(files):

#         # Get info for the files
#         folder   = os.path.dirname(file)
#         basename = os.path.basename(file)
#         basename = os.path.splitext(basename)[0]



#         # Import the file with BioFormats
#         IJ.log("Currently opening " + basename + "...")

#         imps = BFImport(str(file))
#         for imp in imps:

        rm.reset()

        IJ.log("\\Update3:Loading the image from OMERO")
        openImagePlus(HOST, USERNAME, PASSWORD, groupId, image_id)
        # IJ.log("\\Update3:Image loaded")
        imp = IJ.getImage()
        basename = imp.getTitle()


        out_roi        = os.path.join(destination, basename + "_ROI.zip")
        out_single_roi = os.path.join(destination, basename + "_ROI.roi")
        out_3D_rois    = os.path.join(destination, basename + "_3DROIs.zip")
        out_label      = os.path.join(destination, basename + "_nuclei_labels.tif")


        if (os.path.exists(out_roi) or os.path.exists(out_single_roi)):
        # if query_roi.yesPressed():
            IJ.log("\\Update3:Image loaded, ROI found")
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

        nuclei_imp_raw = Duplicator().run(imp,
                                    nuclei_chnl, nuclei_chnl,
                                    1, imp.getNSlices(),
                                    1, imp.getNFrames())

        IJ.log("\\Update3:Processing 1st channel...")
        IJ.run(imp, "Select None", "")

        conA_imp = Duplicator().run(imp,
                                    conA_chnl, conA_chnl,
                                    1, imp.getNSlices(),
                                    1, imp.getNFrames())

        conA_imp_raw = conA_imp.duplicate()


        conA_imp.setTitle("C" + str(conA_chnl))
        IJ.run(conA_imp, "Median...", "radius=2 stack")
        IJ.setAutoThreshold(conA_imp, "Moments dark stack")
        Prefs.blackBackground = True
        IJ.run(conA_imp, "Convert to Mask", "method=Moments background=Dark black")

        IJ.log("\\Update3:Processing 2nd channel...")
        rrna_imp = Duplicator().run(imp,
                                    rrna_chnl, rrna_chnl,
                                    1, imp.getNSlices(),
                                    1, imp.getNFrames())

        rrna_imp_raw = rrna_imp.duplicate()

        rrna_imp.setTitle("C" + str(rrna_chnl))
        IJ.run(rrna_imp, "Median...", "radius=2 stack")
        IJ.setAutoThreshold(rrna_imp, "Moments dark stack")
        Prefs.blackBackground = True
        IJ.run(rrna_imp, "Convert to Mask", "method=Triangle background=Dark black")

        for roi_nb in range(rm.getCount()):
            current_roi = rm.getRoi(roi_nb)
            IJ.log("\\Update3:Currently working on ROI " + current_roi.getName())

            IJ.run(imp, "Select None", "")

            current_label_imp = run_tm(imp, nuclei_chnl, quality_thr, intensity_thr,
                                    circularity_thr, area_thr, cellpose_env, current_roi)
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


        # label_imp.show()
        # sys.exit(0)

        img_nuclei              = ImageInt.wrap(label_imp)
        pop_nuclei              = Objects3DPopulation(img_nuclei)
        nb_nuclei               = pop_nuclei.getNbObjects()
        nucleis_obj_to_remove   = []
        nucleis_index_to_remove = []

        # Loop through objects to remove too big ones
        IJ.log("\\Update3:Removing objects touching borders...")
        for i in range(0, nb_nuclei):
            obj = pop_nuclei.getObject(i)

            # if obj.getVolumeUnit() > max_volume_nuclei:
            #     nucleis_obj_to_remove.append(obj)
            #     nucleis_index_to_remove.append(obj.getValue())
            if obj.touchBorders(img_nuclei, filter_objects_touching_z):
                nucleis_obj_to_remove.append(obj)
                nucleis_index_to_remove.append(obj.getValue())

        for obj in nucleis_obj_to_remove:
            pop_nuclei.removeObject(obj)

        for i in nucleis_index_to_remove:
            IJ.run(label_imp, "Replace/Remove Label(s)",
                        "label(s)=" + str(i) +" final=0")

        # label_imp.show()
        # sys.exit(0)

        dilated_label_imp = label_imp.duplicate()

        nuclei_int           = []
        conA_int             = []
        rrna_int             = []
        dilation_small_array = []

        dilation_steps = max_dilation_radius + 1 - min_dilation_radius

        out_csv = os.path.join(destination,
                        basename + "_results.csv")

        # dilated_label_imp.show()
        # IJ.run(dilated_label_imp, "Dilate labels", "radius=2")
        IJ.log("\\Update3:Dilating the labels...")
        for dilation_radius in range(min_dilation_radius, max_dilation_radius + 1):
            IJ.log("\\Update3: Dilating the labels... Radius " + str(dilation_radius))
            clij2_instance.clear()

            image_title_array.append(basename)
            dilation_array.append(dilation_radius)

            out_c1 = os.path.join(destination,
                                  basename + "_c1_radius_" + str(dilation_radius) + ".tif")
            out_c2 = os.path.join(destination,
                                  basename + "_c2_radius_" + str(dilation_radius) + ".tif")
            out_txt = os.path.join(destination,
                                   basename + "_results_radius_" + str(dilation_radius) + ".txt")


            current_conA_imp = conA_imp.duplicate()
            current_rrna_imp = rrna_imp.duplicate()
            # dilated_label_imp = LabelImages.dilateLabels(dilated_label_imp, dilation_radius)

            dilated_label_imp = dilate_labels_on_gpu(clij2_instance, dilated_label_imp, dilation_radius)

            # dilated_label_imp.show()
            # sys.exit(0)
            dilated_signal_imp = ImageCalculator.run(dilated_label_imp, label_imp,
                            "Subtract create stack")

            # test.show()
            # sys.exit(0)

            # ImageCalculator.run(dilated_label_imp, label_imp,
            #                 "Subtract stack")

            dilated_conA_imp = ImageCalculator.run(dilated_signal_imp, current_conA_imp,
                                "AND create stack")
            dilated_conA_imp.setTitle("Dilated_conA")
            dilated_rrna_imp = ImageCalculator.run(dilated_signal_imp, current_rrna_imp,
                                "AND create stack")
            dilated_rrna_imp.setTitle("Dilated_rrna")

            # current_conA_imp.show()
            # current_rrna_imp.show()
            # dilated_conA_imp.show()
            # dilated_rrna_imp.show()

            dilated_conA_imp_thresh = dilated_conA_imp.duplicate()
            dilated_rrna_imp_thresh = dilated_rrna_imp.duplicate()

            IJ.setRawThreshold(dilated_conA_imp_thresh, 1, 65535)
            IJ.setRawThreshold(dilated_rrna_imp_thresh, 1, 65535)

            IJ.run(dilated_conA_imp_thresh, "Convert to Mask",
                   "method=Default background=Dark black")

            IJ.run(dilated_rrna_imp_thresh, "Convert to Mask",
                   "method=Default background=Dark black")

            dilated_conA_imp_thresh.show()
            dilated_rrna_imp_thresh.show()
            # label_imp.show()

            # sys.exit(0)


            IJ.log("\\Update3:Running JACoP...")
            IJ.run("JACoP ", "imga=" + dilated_conA_imp_thresh.getTitle() +
                " imgb=" + dilated_rrna_imp_thresh.getTitle() +" thra=1 thrb=1 pearson mm")

            dilated_conA_imp_thresh.close()
            dilated_rrna_imp_thresh.close()
            current_conA_imp.close()
            current_rrna_imp.close()



            IJ.log("\\Update3:Saving images...")
            if save_tiffs:
                BFExport(label_imp, out_label)
                BFExport(dilated_conA_imp,out_c1)
                BFExport(dilated_rrna_imp, out_c2)

            img_nuclei = ImageInt.wrap(label_imp)
            img_conA   = ImageInt.wrap(dilated_conA_imp)
            img_rrna   = ImageInt.wrap(dilated_rrna_imp)

            pop_nuclei = Objects3DPopulation(img_nuclei)
            pop_conA   = Objects3DPopulation(img_conA)
            pop_rrna   = Objects3DPopulation(img_rrna)

            nb_nuclei  = pop_nuclei.getNbObjects()



            ih_nuclei = ImageHandler.wrap(nuclei_imp_raw)
            ih_conA   = ImageHandler.wrap(conA_imp_raw)
            ih_rrna   = ImageHandler.wrap(rrna_imp_raw)

            for i in range(0, nb_nuclei):
                current_nuc  = pop_nuclei.getObject(i)
                nuclei_int.append(current_nuc.getPixMeanValue(ih_nuclei))

                nuclei_value = current_nuc.getValue()

                current_conA = pop_conA.getObjectByValue(nuclei_value)
                if current_conA:
                    conA_int.append(current_conA.getPixMeanValue(ih_conA))
                else:
                    conA_int.append(0)

                current_rrna = pop_rrna.getObjectByValue(nuclei_value)
                if current_rrna:
                    rrna_int.append(current_rrna.getPixMeanValue(ih_rrna))
                else:
                    rrna_int.append(0)

            dilated_conA_imp.close()
            dilated_rrna_imp.close()
            dilated_signal_imp.close()
            dilated_label_imp.close()

            dilation_small_array.extend([dilation_radius] * nb_nuclei)


            # IJ.selectWindow("Log")
            # IJ.saveAs("Text", out_txt)
            log_text = IJ.getLog()
            log_array = log_text.split("\n")

            pearson_index = [ i for i, line in enumerate(log_array) if 'Pearson' in line ]
            manders_index = [ i for i, line in enumerate(log_array) if 'Manders' in line ]

            pearson_array.append(log_array[pearson_index[0] + 1][2:])
            m1_array.append(log_array[manders_index[0] + 1][3:8])
            m2_array.append(log_array[manders_index[0] + 2][3:8])

            IJ.log("\\Clear")
            IJ.run("Collect Garbage", "")

        with open(out_csv, 'wb') as f:
            writer = csv.writer(f, delimiter = ';')
            writer.writerow(
                ["Image name; Dilation, Nuclei ID; Nuclei intensity; C1 intensity; C2 intensity"])
            writer.writerows(izip([basename] * nb_nuclei * dilation_steps,
                                  dilation_small_array,
                                  range(1, nb_nuclei + 1) * dilation_steps,
                                  nuclei_int, conA_int, rrna_int))


        imp.close()
        label_imp.close()
        nuclei_imp_raw.close()
        conA_imp.close()
        conA_imp_raw.close()
        rrna_imp.close()
        rrna_imp_raw.close()






    with open(coexpr_csv, 'wb') as f:
        writer = csv.writer(f, delimiter =';')
        writer.writerow(
            ["Image Name; Dilation radius; Pearson coefficient; M1 coefficient; M2 coefficient"]
        )
        writer.writerows(izip(image_title_array, dilation_array, pearson_array,
                              m1_array, m2_array))

	IJ.log("\\Clear")
    IJ.log("\\Update2:Script DONE ! ")
    IJ.log("\\Update3:")

finally:
    gateway.disconnect()


