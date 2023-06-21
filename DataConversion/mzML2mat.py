from DataConvert import DataConvert
import time

RawDataFile = '20220624-004-MBSI-neg-ms1-slic-im3-tran40.mzML'
FileName = "PixelData"
start = time.clock()
DataConvert(RawDataFile,FileName,200)
elapsed = (time.clock() - start)
print("Time used:",elapsed)