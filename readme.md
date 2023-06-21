# 1. Introduction

Mass spectrometry imaging (MSI) is widely used for in situ ionization and multiplexed, label-free detection of metabolites in cells and tissues. One of the key challenges in MSI is the structural annotation of the hundreds of detected peaks. We have developed a spatial structural metabolomics method enabled by ion mobility–data independent acquisition coupled with MSI (IM-DIA imaging). This set of codes are used for data processing, MS/MS deconvolution, and image reconstruction of data produced by IM-DIA imaging.

The workflow of IM-DIA imaging is shown in **Fig. 1a**.The IM-DIA acquisition method consists of two IM-MS/MS functions using different post-IM collision energy (Fig. 1a), namely HDMS<sup>E</sup> in MassLynx. The IM-MS data acquired using these two functions are named as low energy (LE) spectra and high energy (HE) spectra, respectively. LE and HE spectra were acquired on adjacent pixels across the whole imaging area, yielding a multi-dimensional MS data indexed by pixel position, drift arrival time, and m/z.

![IM-DIA imaging workflow and data processing](/Figures/Fig%201.png)
<small>**Fig. 1 Scheme of the MMSD imaging data acquisition workflow and MS/MS deconvolution.** **a** Ion flow on the Cyclic-IMS. Ions produced from a pixel are separated by the cyclic IMS cell, followed by co-fragmentation in the collision cell and mass detection in the TOF. **b** IM-MS profiles of intact precursors (LE spectra) and IM-MS/MS profiles of fragments (HE spectra) were acquired on adjacent pixels using low and high collision energies, yielding IM-MS spectra and IM-MS/MS spectra, respectively. **c~f** Workflow of the MS/MS deconvolution algorithm. Drift time distribution profiles of precursors (P) and fragments (F) were extracted from LE and HE spectra, respectively (**c**). The fragment ion intensity is assumed to be the linear combination of the contributions from all precursors ($F_j=∑_iP_i⋅c_{ij}$). The fragmentation coefficients (c) were solved by nonnegative least square optimization ($\overline{C}=\arg_{C>=0}min$), yielding a fragmentation tree containing fragmentation coefficients from each precursor to each fragment (**d**). The MS/MS deconvolution can be either applied to ROI averaged data, yielding pseudo MS/MS spectra for molecule identification (**e**), or performed on each pixel for multiplexed MS/MS imaging (**f**).</small>

For data processing and spectrum deconvolution, the Raw MS data files were first converted to MATLAB data format (Methods). Peak picking of precursors and fragments was performed on averaged LE and HE MS spectra. Then ion drift time profiles were extracted for each detected precursor and fragment (**Fig. 1c**). It is reasonable to assume that the conversion ratio of precursor ions to different fragments (named fragmentation coefficient, *c*) is constant under a constant collision energy, independent of the precursor ion intensity. The deconvolution algorithm fits the drift time profiles of each fragment ion as the linear combination of drift time profiles of a set of precursors by non-negative least square optimization, yielding a fragmentation tree showing the fragmentation coefficients of different precursors to different fragments (**Fig. 1d**). Conversely, a pseudo-MS/MS spectrum can be generated for each precursor. The MS/MS deconvolution can either be performed on the averaged data across a selected region of interest (ROI) (ROI deconvolution) for accurate molecular identification (**Fig. 1e**), or performed on all pixels on a tissue section (pixel deconvolution) for untargeted structural MSI via multiplexed MS/MS imaging (**Fig .1f**).

# 2. Quick start

## 2.1 Requirements

- Software: **MassLynx V4.2**, **HDIimaging V1.5**, **Python 3.7 (pymzml, scipy, numpy)**, **MATLAB R2021b**
- Hardware: All codes were successfully run on a laptop (Windows 10) with a **11th Gen Intel(R) Core(TM) i5-11400H @ 2.70GHz CPU** and **16 GB RAM**.

## 2.2 Data conversion

Data conversion takes most of the time. Data compression is recommended if the raw data file is too large.

- Waters .raw data file to .mzML data file : **MSConvert**

- .mzML to .mat : run **\DataConversion\mzML2mat.py**

    <small>**Note**: Change the RawDataFile and DTScanNumber (defalt=200) if necessary. </small>

Output .mat files are stored in **\DataConversion\Output_mat_DataFile**

## 2.3 Data pre-processing

Perform MS re-calibration if necessary.

## 2.4 Data overview using HDImaging

- Generate .txt files of the imaging data using **HDIimaging**
- Determine the number of rows and columns of LE and HE pixels
- Generate .txt files of the averaged mass spectra of LE and HE data using **MassLynx**

## 2.5 Spatial segmentation

- Run **\Spatial Segmentation\Segmentation.m**
  
  ```matlab
    % load("Mask.mat");
     Mask = ones(R,C1+C2);
  ```

  First, perform a pre-segmentation to determine the tissue area, generate **Mask.mat**. Background pixels are 0 and tissue pixels are 1.

- Run **\Spatial Segmentation\Segmentation.m**

  ```matlab
    load("Mask.mat");
    % Mask = ones(R,C1+C2);
  ```

  Segmentation result is stored in **ClassID.mat**

## 2.6 ROI deconvolution

- Run **\Deconvolution\SegDataAverage.m**

    ```matlab
    ParentIonList = ILG.Relative(RawParent,700,900,0.01);
    FragmentIonList = ILG.Relative(RawFragment,100,900,0.005);
    % ParentIonList = importdata("PixelDecFeature\ParentIonListMSI.txt");
    % FragmentIonList = importdata("PixelDecFeature\FragmentIonListMSI.txt");
    ```

    <small>**Note**: Change file path, file name, parameters if necessary.</small>

## 2.7 Lipid annotation

- MS1 search: Copy **ParentIonList.mat** to **\Annotation** and run **MS1Searchtst.m**

<small>**Note**: Change calibration parameters as per your experiment.</small>

- MS2 search: Copy **ParentIonList.mat**, **FragmentIonList.mat**, and **CombinedFragmentationMap.mat** to **\Annotation** and run **MS2Searchtst.m**

## 2.8 Pixel deconvolution

- Determine annotated parent ion list and fragment ion list in **2.7**

- (Optional) Run **\Deconvolution\SegDataAverage.m**

    ```matlab
    %ParentIonList = ILG.Relative(RawParent,700,900,0.01);
    %FragmentIonList = ILG.Relative(RawFragment,100,900,0.005);
    ParentIonList = importdata("PixelDecFeature\ParentIonListMSI.txt");
    FragmentIonList = importdata("PixelDecFeature\FragmentIonListMSI.txt");
    ```

    <small>**Note**: Change file path, file name, parameters if necessary. This step is optional. But reduce the number of targeted precursors and fragments should largely save the run time.</small>

- Choose Parent-Fragment pair for MS/MS imaging (**FragmentationFeature.txt**)

- Run **\Deconvolution\MSIResconstruction2.m**
  
    <small>**Note**: Change file path, file name, parameters if necessary.</small>