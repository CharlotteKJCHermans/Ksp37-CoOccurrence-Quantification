# ============================================================
# CO-OCCURRENCE of Ksp37 and subcellular markers (protein of 
# interest, POI)
# - Segments cells on channel 0
# - Filters to Ksp37+ and POI+ double-positive cells only
# (channel mean signal > threshold)
# - Saves a table with the co-occurrence results per cell
# - Saves a JPEG outline of double positive cells
# ============================================================

from ij import IJ, WindowManager, Prefs
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.macro import Interpreter
from sc.fiji.coloc import Coloc_2
from java.io import File
from java.awt import Color
from ij.gui import Overlay
import java.lang.System as JSystem
import os

# --- Helpers ---
def get_field(obj, name):
    cls = obj.getClass()
    c = cls
    f = None
    while c is not None:
        try:
            f = c.getDeclaredField(name)
            break
        except:
            c = c.getSuperclass()
    if f is None:
        raise Exception("Field not found: " + name)
    f.setAccessible(True)
    return f.get(obj)

def get_value(results, name):
    for vr in results.values():
        if vr.name == name:
            if vr.isNumber:
                return vr.number
            else:
                return vr.value
    return float('nan')

# ============================================================
# 0. Preparation
# ============================================================
inputDir = IJ.getDirectory("Choose the folder containing your .nd2 files")
mainOutputDir = IJ.getDirectory("Choose your output folder")

subFolder = mainOutputDir + "Cell_segmentation" + File.separator
subDir = File(subFolder)
if not subDir.exists():
    subDir.mkdirs()

fullTablePath = mainOutputDir + "Results_CoOccurrence.csv"

tableFile = File(fullTablePath)
if tableFile.exists():
    rt = ResultsTable.open(fullTablePath)
else:
    rt = ResultsTable()

ksp37_radius = 6 
POI_radius = 7
signalThreshold_Ksp37 = 8 
signalThreshold_POI = 20
ksp37_min_offset = 30
ksp37_max_offset = 180
POI_min_offset = 50
POI_max_offset = 280

# ============================================================
# Per-image processing
# ============================================================
def process_image(nd2Path):
    IJ.log("Opening: " + nd2Path)
    IJ.run("Bio-Formats Importer",
           "open=[" + nd2Path + "] color_mode=Default view=Hyperstack "
           "stack_order=XYCZT split_channels")

    title0 = None
    title1 = None
    title2 = None
    for t in WindowManager.getImageTitles():
        if t.endswith("C=0"):
            title0 = t
        if t.endswith("C=1"):
            title1 = t
        if t.endswith("C=2"):
            title2 = t

    if title0 is None or title1 is None or title2 is None:
        IJ.log("WARNING: could not find all 3 channels for " + nd2Path + " - skipping")
        IJ.run("Close All")
        return

    imp0 = WindowManager.getImage(title0)
    imp1 = WindowManager.getImage(title1)
    imp2 = WindowManager.getImage(title2)

    rm = RoiManager.getInstance()
    if rm is None:
        rm = RoiManager()
    rm.reset()

    # 1. Subtract Background
    IJ.run(imp0, "Subtract Background...", "rolling=60 separate")
    IJ.run(imp1, "Subtract Background...", "rolling=" + str(ksp37_radius) + " separate")
    IJ.run(imp2, "Subtract Background...", "rolling=" + str(POI_radius) + " separate")

    # 2. Create mask of membrane stain
    maskImp = imp0.duplicate()
    maskImp.setTitle("maskCells")
    IJ.run(maskImp, "Gaussian Blur...", "sigma=2")
    IJ.setAutoThreshold(maskImp, "Huang dark")
    Prefs.blackBackground = True
    IJ.run(maskImp, "Convert to Mask", "")

    # 3. Analyze particles to create ROIs
    IJ.run(maskImp, "Analyze Particles...", "size=20-Infinity include add")
    IJ.setForegroundColor(255, 255, 255)
    WindowManager.setTempCurrentImage(maskImp)
    rm.runCommand("Fill")
    WindowManager.setTempCurrentImage(None)

    nAll = rm.getCount()
    keepIndices = []
    for i in range(nAll):
        imp1.setRoi(rm.getRoi(i))
        mean_Ksp37 = imp1.getStatistics().mean
        imp2.setRoi(rm.getRoi(i))
        mean_POI = imp2.getStatistics().mean
        if mean_Ksp37 > signalThreshold_Ksp37 and mean_POI > signalThreshold_POI:
            keepIndices.append(i)

    allRois = rm.getRoisAsArray()
    rm.reset()
    for idx in keepIndices:
        rm.addRoi(allRois[idx])

    imp1.killRoi()
    imp2.killRoi()

    if rm.getCount() == 0:
        IJ.log("No double-positive cells found for " + nd2Path + " - skipping")
        IJ.run("Close All")
        return

    # 4. Overlay image
    IJ.run(imp1, "Grays", "")
    IJ.run(imp2, "Grays", "")
    IJ.setMinAndMax(imp1, imp1.getProcessor().getMin() + ksp37_min_offset,
                          imp1.getProcessor().getMin() + ksp37_max_offset)
    IJ.setMinAndMax(imp2, imp2.getProcessor().getMin() + POI_min_offset,
                          imp2.getProcessor().getMin() + POI_max_offset)

    IJ.run("Merge Channels...", "c2=[" + title1 + "] c6=[" + title2 + "] create keep ignore")
    overlayImp = IJ.getImage()
    overlayImp.setTitle("Coloc2_Overlay")

    ov = Overlay()
    for roi in rm.getRoisAsArray():
        roi.setStrokeColor(Color.WHITE)
        roi.setStrokeWidth(1)
        ov.add(roi)
    overlayImp.setOverlay(ov)

    flatImp = overlayImp.flatten()
    IJ.run(flatImp, "Scale Bar...", "width=20 height=10 font=32 horizontal bold")
    IJ.saveAs(flatImp, "Jpeg", subFolder + title1 + "_Coloc2_overlay.jpg")
    overlayImp.close()
    flatImp.close()

    # 5. Measure colocalization per cell
    rm.deselect()
    c2 = Coloc_2()
    ok = c2.initializeSettings(imp1, imp2, 3, 0, False, False, False, False, False, False,
                                True, True, False, False, True, 2, 100)
    if not ok:
        IJ.log("initializeSettings() failed for " + nd2Path + " - skipping")
        IJ.run("Close All")
        return

    masks = get_field(c2, "masks")
    img1 = get_field(c2, "img1")
    img2 = get_field(c2, "img2")

    JSystem.setProperty("java.awt.headless", "true")

    for i in range(len(masks)):
        IJ.log("Cell " + str(i + 1) + " / " + str(len(masks)) + " done")
        mi = masks[i]
        roi = get_field(mi, "roi")
        mask = get_field(mi, "mask")
        results = c2.colocalise(img1, img2, roi, mask, None)

        currentRow = rt.getCounter()
        rt.incrementCounter()
        rt.setValue("Image_Name", currentRow, title1)
        rt.setValue("Cell ID", currentRow, i + 1)
        rt.setValue("Ksp37_BG_Radius", currentRow, ksp37_radius)
        rt.setValue("POI_BG_Radius", currentRow, POI_radius)
        rt.setValue("Signal Threshold Ksp37", currentRow, signalThreshold_Ksp37)
        rt.setValue("Signal Threshold POI", currentRow, signalThreshold_POI)
        rt.setValue("Ksp37_MinOffset", currentRow, ksp37_min_offset)
        rt.setValue("Ksp37_MaxOffset", currentRow, ksp37_max_offset)
        rt.setValue("POI_MinOffset", currentRow, POI_min_offset)
        rt.setValue("POI_MaxOffset", currentRow, POI_max_offset)
        rt.setValue("Coloc Threshold Ksp37", currentRow, get_value(results, "Ch1 Max Threshold"))
        rt.setValue("Coloc Threshold POI", currentRow, get_value(results, "Ch2 Max Threshold"))
        rt.setValue("Costes P-Value", currentRow, get_value(results, "Costes P-Value"))
        rt.setValue("Pearson R", currentRow, get_value(results, "Pearson's R value (no threshold)"))
        rt.setValue("Pearson R (above threshold)", currentRow, get_value(results, "Pearson's R value (above threshold)"))
        rt.setValue("Manders M1 (above zero)", currentRow, get_value(results, "Manders' M1 (Above zero intensity of Ch2)"))
        rt.setValue("Manders M2 (above zero)", currentRow, get_value(results, "Manders' M2 (Above zero intensity of Ch1)"))
        rt.setValue("Manders M1 (autothreshold)", currentRow, get_value(results, "Manders' tM1 (Above autothreshold of Ch2)"))
        rt.setValue("Manders M2 (autothreshold)", currentRow, get_value(results, "Manders' tM2 (Above autothreshold of Ch1)"))
        rt.setValue("Spearman", currentRow, get_value(results, "Spearman's rank correlation value"))

    JSystem.setProperty("java.awt.headless", "false")

    # Save after every image - protects progress if a later file fails
    rt.save(fullTablePath)
    IJ.run("Close All")

# ============================================================
# Batch loop over all .nd2 files in the input folder
# ============================================================
Interpreter.batchMode = True

nd2Files = [f for f in os.listdir(inputDir) if f.endswith(".nd2")]
nd2Files.sort()
IJ.log("Found " + str(len(nd2Files)) + " .nd2 files")

for fname in nd2Files:
    process_image(inputDir + fname)

Interpreter.batchMode = False
rt.show("Results_CoOccurrence")
IJ.run("Clear Results")