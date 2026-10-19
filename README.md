# Ksp37-CoOccurrence-Quantification
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23474924.svg)](https://doi.org/10.5281/zenodo.23474924)

Fiji (Jython) script for quantifying co-occurrence of Ksp37 with a protein of interest (POI) in single cells, for example a vesicle marker or an effector molecule. The script can also be used for other proteins, provided the proteins to analyze are in channel 1 and channel 2 and the cell mask is in channel 0 (see User parameters).

## Description
Cells are segmented on a membrane stain. The mean Ksp37 and POI signal of each cell is compared with manually set thresholds to identify double-positive cells (expressing both Ksp37 and POI). Co-occurrence statistics are then calculated for these cells with Coloc 2. The script saves a results table and an outline image (JPEG) highlighting the double-positive cells included in the analysis.

## Requirements
- Fiji (ImageJ 1.54p; Java 21.0.7, 64-bit), macOS Tahoe 26.6.2
- Bio-Formats plugin (version 8.5.0), used to open .nd2 files
- Coloc 2 plugin (version 3.1.0), used for colocalization statistics
- Both plugins are included in the default Fiji installation. No additional update sites required.
- The script is written in Python (Jython) and runs in Fiji's Script Editor.

## Input format
- .nd2 files
- Channel 0: cell mask channel used for segmentation
- Channel 1: Ksp37 staining
- Channel 2: subcellular marker (POI)
- Files without all three channels are skipped.

## Example data
Sample images (.nd2) are available on Zenodo: https://doi.org/10.5281/zenodo.23473247. 
Download them into a folder and run the script with that folder as input. Compare your results with `Ksp37-CoOccurrence-Quantification-ExampleData/Expected_Output`.

## How to run
1. In Fiji, click File → New → Script…
2. In the Script Editor, click File → Open… and select the .py file (or drag the file onto Fiji).
3. Check that Language → Python is selected.
4. Click Run, then choose the input folder with the images and the output folder.

The results table is saved in the output folder, and a subfolder `Cell_segmentation` is created for the overlay images.

## Segmentation and analysis (fixed in the script)
- Background subtraction, channel 0: rolling ball radius 60
- Gaussian blur (sigma 2), then Huang auto-threshold and Analyze Particles with a minimum size of 20
- Scale bar: 20 µm, requires calibrated images
- Coloc 2 settings: Channel 1, Channel 2, ROI manager, Costes threshold regression, Spearman's Rank Correlation, Manders' Correlation, Costes' Significance Test, PSF of 2.0, 100 Costes randomizations

## User parameters
Set in the script before the per-image processing.
- `ksp37_radius`, `POI_radius`: rolling ball radii for background subtraction of channels 1 and 2. Values are optimized to remove cytoplasmic background while avoiding loss of true signal.
- `signalThreshold_Ksp37`, `signalThreshold_POI`: mean signal (after background subtraction) above which a cell counts as positive for that channel. Only cells above both thresholds are analyzed.
- `ksp37_min_offset`, `ksp37_max_offset`, `POI_min_offset`, `POI_max_offset`: display range of the overlay image, relative to the minimum pixel value of each channel. They do not affect the quantification.

Values were tuned for our stainings. Adjust all of them for your own data. To quantify other proteins, put them in channel 1 and channel 2 and adjust `ksp37_radius`, `POI_radius`, `signalThreshold_Ksp37`, `signalThreshold_POI` and the display offsets to your signal. Note that output column and file names keep the labels “Ksp37” and “POI”.

## Outputs
- `Results_CoOccurrence.csv`: one row per double-positive cell. Re-running on the same output folder adds new rows and does not overwrite earlier ones, so use an empty output folder for each analysis to avoid duplicates.
- `Cell_segmentation/<image name>_Coloc2_overlay.jpg`: overlay of the Ksp37 and marker channels with the analyzed cells outlined.

CSV columns:
- `Image_Name`, `Cell ID`
- `Ksp37_BG_Radius`, `POI_BG_Radius`, `Signal Threshold Ksp37`, `Signal Threshold POI`, `Ksp37_MinOffset`, `Ksp37_MaxOffset`, `POI_MinOffset`, `POI_MaxOffset`: parameter values used
- `Coloc Threshold Ksp37`, `Coloc Threshold POI`: thresholds determined by Coloc 2 for each channel
- `Costes P-Value`
- `Pearson R`, `Pearson R (above threshold)`: Pearson's R value without and with Coloc 2 threshold applied, respectively
- `Manders M1/M2 (above zero)`, `Manders M1/M2 (autothreshold)`: M1 and M2 as defined by Coloc 2 (see its documentation)
- `Spearman`: Spearman's rank correlation

## License
MIT

## Cite
If you use this script, please cite:

Hermans, CKJC (2026). Ksp37-CoOccurrence-Quantification: Fiji script for co-occurrence of Ksp37 with a protein of interest (v1.0.0).
Zenodo. https://doi.org/10.5281/zenodo.23474924

Paper: [link to be added once published]
