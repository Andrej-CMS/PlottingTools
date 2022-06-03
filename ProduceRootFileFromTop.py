import os
import sys
import ROOT
import ROOT.TH1D as TH1D
import ROOT.TH2D as TH2D
import ROOT.TGraph2D as TGraph2D
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TF1, TString
from ROOT import kBlack, kGreen, kRed
from ROOT import TFile, gStyle, TRatioPlot, gPad, TGraphErrors, TGraphAsymmErrors, TObject
import math
import array
import glob
import subprocess
import os
import numpy as np
import pandas as pd

from optparse import OptionParser
from plottingUtils import *


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename", 
                    help="Input file name", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile",
                    help="Output file name", metavar="FILE")

    (options, args) = parser.parse_args()
    
    outputPath = "HistoListsFromTextFiles"
    outputFileName = options.outputFile
    debug = False
    
    fileNameList = options.filename.split(",")
    
    binCenters = []
    binContent = [] 
    binUncertainties = []
    histoList = []
    
    fileContentDictionary = {}

    #read input files
    #all input files should be given
    #the script creates all histograms in one single file
    #this way, all variations can be read from one output ROOT file
    
    #read in all lines from input top files
    #structure is: [bin-center, bin-content, bin-uncertainties]
    #fill content into dictionary
    for fileName in fileNameList:
        with open(fileName) as file:

            del binCenters[:]
            del binContent[:]
            del binUncertainties[:]

            for line in file:
                lineContent = line.split()
                
                binCenters.append(float(lineContent[0]))
                binContent.append(float(lineContent[1]))
                binUncertainties.append(float(lineContent[2]))
                
            fileContentDictionary.update({fileName : [binCenters, binContent, binUncertainties]})
                
    
    if debug:
        print(fileContentDictionary)
        print(binCenters)
        print(binContent)
        print(binUncertainties)
    
    #contruct ROOT TH1 histograms from input dictionary
    #Name of the histogram: mass_PDF_scale_distribution
    for key in fileContentDictionary:
        fileName = makeFileNameFromPath(key)
        
        histoList.append(makeHistoFromLists(fileContentDictionary.get(key)[0], 
                                            fileContentDictionary.get(key)[1],
                                            fileContentDictionary.get(key)[2], 
                                            fileName))
        
    #checked manually that the numbers are exactly the same as in the input files
    if debug:
        for xBin in range(histoList[0].GetNbinsX()):
            print(xBin, histoList[0].GetBinCenter(xBin))
          
    #create output folder and write output file
    outputList = TList()
    for histo in histoList:
        outputList.Add(histo)
        
    if not os.path.exists(outputPath):
        # Create a new directory because it does not exist 
        os.makedirs(outputPath)
        print("The new directory is created!")
        
    outputFile = TFile(outputPath+"/"+outputFileName+".root","RECREATE");
    outputList.Write("histoList", TObject.kSingleKey);
    outputFile.Close()
    
    print("Finished processing ", outputFileName)
    