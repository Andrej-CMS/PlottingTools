import os
import sys
import ROOT
import ROOT.TH1D as TH1D
import ROOT.TH2D as TH2D
import ROOT.TGraph2D as TGraph2D
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TF1, TString
from ROOT import kBlack, kGreen, kRed
from ROOT import TFile, gStyle, TRatioPlot, gPad, TGraphErrors, TGraphAsymmErrors
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
    
    #read input files
    #all input files should be given
    #the script creates all histograms in one single file
    #this way, all variations can be read from one output ROOT file
    
    #read in all lines from input top files
    #structure is: [bin-center, bin-content, bin-uncertainties]
    #fill content into dictionary
    
    #contruct ROOT TH1 histograms from input dictionary
    #Name of the histogram: mass_PDF_scale_distribution
          
    #create output folder and write output file
    
    