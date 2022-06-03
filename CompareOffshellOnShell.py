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
    parser.add_option("-f", "--onShellFile", dest="onShellFile", 
                    help="Input file name", metavar="FILE")
    parser.add_option("-g", "--offShellFile", dest="offShellFile",
                    help="Path to input file with ", metavar="FILE")

    (options, args) = parser.parse_args()
    
    outputPath  = "OnShell_Offshell_comparison"

    numberOfBins = 5
    xBinsLow    = 0.0
    xBinsHigh    = 1.0
    xBins = array.array('d',np.linspace(xBinsLow, xBinsHigh, numberOfBins+1))
    # get files

    onShellFile     = ROOT.TFile.Open( options.onShellFile ," READ ")
    offShellFile    = ROOT.TFile.Open( options.offShellFile ," READ ")

    # Make histograms
    # For Offshell and Onshell each three histograms

    variable    = "27_rho"
    PDFscale    = "NNPDF3_HT_6"
    mass       = "mt173"
    variableOnShell = "rho_fb(ttbar+jet)>_0"
    

    offShellHistoList   = offShellFile.Get("histoList")
    offShellHisto       = offShellHistoList.FindObject(variable+"_"+PDFscale+"_"+mass)
    print("offShellHisto Integral ", offShellHisto.Integral())
    onShellHisto        = onShellFile.Get(variableOnShell)



    offShellHisto   = offShellHisto.Rebin(numberOfBins, "hnew offshell", xBins)
    onShellHisto    = onShellHisto.Rebin(numberOfBins, "hnew onshell", xBins)
    # Plot histograms on one canvas with ratio
    canvas, pad1, pad2 = createCanvasPads()
    canvas.cd()
    pad1.cd()
    
    #create legend
    legend = ROOT.TLegend(.58,.65,.77,0.87)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)
    
    offShellHisto.Scale(1./offShellHisto.Integral())
    offShellHisto.Draw()
    
        
    onShellHisto.SetLineColor(kRed)
    onShellHisto.Scale(1./onShellHisto.Integral())
    onShellHisto.Draw("same")
    
    pad2.cd()
    
    #we want to calculate DeltaR/R * 1/DeltaM

    ratio = createRatio(offShellHisto, onShellHisto, "Top quark p_{T}", "#frac{Pythia 8}{Herwig 7}")
    
    ratio.Draw()
    # Create a new directory because it does not exist 
    
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
        print("The new directory is created!")

    canvas.Update()
    string_normalize = ""
    # if normalize:
    #     string_normalize = "_normalized"

    canvas.SaveAs(outputPath+"/"+variable+"_"+PDFscale+"_"+mass+".pdf")
    