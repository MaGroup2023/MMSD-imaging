
def DataConvert(RawDataFile,FileName,DTScanNumber):
    import pymzml
    import numpy as np
    from scipy.io import savemat
    
    RawData = pymzml.run.Reader(RawDataFile)
    PixelData = [[]]*DTScanNumber #Data in a pixel

    for i, spec in enumerate(RawData):
        FunctionNumber = spec.id_dict['function']
        DTIndex = (spec.id_dict['scan'])%DTScanNumber #Current DT index
        CurrentPixel = i//(2*DTScanNumber)
        PixelData[DTIndex-1] = spec.peaks("raw")
        if DTIndex == 0:
            if FunctionNumber == 1: #Low energy
                path = "Output_mat_DataFile\\LE"+FileName+str(CurrentPixel+1).zfill(6)+".mat"
                savemat(path,{'PixelData':PixelData})
            else: #High energy
                path = "Output_mat_DataFile\\HE"+FileName+str(CurrentPixel+1).zfill(6)+".mat"
                savemat(path,{'PixelData':PixelData})



