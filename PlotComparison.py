from ast import For
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

from plottingUtils import *

ROOT.gROOT.SetBatch(True)


outputPath = "ComparisonPlots"
variable = "weight_mc"
normalize = False


# Set properties of plots based on dictionary of variables
variableFixedOrder = getListOfVariableProperties(variable)[0]
xTitle = getListOfVariableProperties(variable)[1]
yTitle = getListOfVariableProperties(variable)[2]
numberOfBins = getListOfVariableProperties(variable)[3]
xBinsLow = getListOfVariableProperties(variable)[4][0]
xBinsHigh = getListOfVariableProperties(variable)[4][1]
xBins = array.array('d',np.linspace(xBinsLow, xBinsHigh, numberOfBins+1))

print xBins

ratioMin = getListOfVariableProperties(variable)[5][0]
ratioMax = getListOfVariableProperties(variable)[5][1]
histoTitle = getListOfVariableProperties(variable)[6]
totalXS_FixedOrder_Name = getListOfVariableProperties(variable)[7]
treeName = getListOfVariableProperties(variable)[8]
    

yTitleRatio = "ratio"
yTitleAddition = ""
if normalize:
    if variable == "weight_mc":
        yTitleAddition = "Normalized "
    else:
        yTitleAddition = "#frac{1}{#sigma}"

histoList = []
histo     = ROOT.TH1D()
legendNameList = []


# run over all input files and plot them into one histogram
for argNumber in range(len(sys.argv)-1):
    
    fileName = sys.argv[argNumber+1]
    print " Reading from file", fileName 
    inFile = ROOT.TFile.Open( fileName ," READ ")

    legendNameList.append(fileName[:-5])

    tree  = inFile.Get (treeName)

    histo = ROOT.TH1D( fileName, histoTitle+str(argNumber+1) , numberOfBins , xBins )
    # SetDirectory(0) to keep the histogram instead of being deleted after file close
    histo.SetDirectory(0)
    histo.SetMinimum(0.001)

    sumOfWeights = 0.
    factor_GeV = 1.
    isRhoVariable = True
    
    if (not "_rho_" in variable):
        factor_GeV = 1000.
        isRhoVariable = False
    elif not "weight" in variable:
        factor_GeV = 1.

    for entryNum in range (0 , tree.GetEntries()):

        tree.GetEntry( entryNum )

        if variable == "jet_HT":
            eventValue      = getattr( tree ,"jet_pt")
        else:
            eventValue      = getattr( tree ,variable)

        eventWeight     = getattr( tree , "weight_mc")
        sumOfWeights    += eventWeight
        jet_HT = 0.

        # if isRhoVariable:
        #     jetPt = getattr( tree ,"MC_j_afterFSR_pt")
        #     if jetPt < 30.:
        #         continue
        if variable == "jet_pt":
            for value in eventValue:
                histo.Fill(value/factor_GeV, eventWeight)
                
        elif variable == "jet_HT":
            for value in eventValue:
                jet_HT+=value
            histo.Fill(jet_HT/factor_GeV, eventWeight)
        elif variable == "weight_mc": 
            histo.Fill(eventValue)
        else:
            histo.Fill(eventValue/factor_GeV, eventWeight)

    crossSectionFactor = sumOfWeights/tree.GetEntries()

    if normalize:
        histo.Scale(1./histo.Integral())
    elif not "weight" in variable:
        #need to normalize to expected events at 1 inverse femtobarn
        histo.Scale(1./crossSectionFactor)

    print "crossSectionFactor, histo   ", crossSectionFactor, " ", histo

    histoList.append(histo)


#print histoList

# set style for canvas
gStyle.SetPadLeftMargin(0.16)
gStyle.SetTitleFontSize(0.05)

canvas, pad1, pad2 = createCanvasPads()
canvas.cd()

#create legend
legend = ROOT.TLegend(.58,.65,.77,0.87)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.05)

pad1.cd()
gPad.SetLogy(0)

listOfLegendNames = []
# loop over all histos in histolist and draw them into the canvas and ratios
for histoNumber in range(len(histoList)):
    
    histoList[histoNumber].GetYaxis().SetTitle(yTitleAddition+yTitle)
    histoList[histoNumber].GetYaxis().SetTitleSize(0.05)
    histoList[histoNumber].GetXaxis().SetTitle(xTitle)
    histoList[histoNumber].SetMaximum(1.3*histoList[argNumber].GetMaximum())
    histoList[histoNumber].SetTitle(histoTitle)

    histoList[histoNumber].SetStats(0)
    histoList[histoNumber].SetLineWidth(2)
    histoList[histoNumber].SetLineColor(histoNumber+1)

    legendName_tmp  = legendNameList[histoNumber].split("/")
    legendName      = legendName_tmp[len(legendName_tmp)-1]

    legend.AddEntry(histoList[histoNumber], getLegendNames(legendName),"L")
    listOfLegendNames.append(legendName)
    
    #calculate fraction of negative weights
    if variable == "weight_mc":
        print "Fraction of negative events ", getLegendNames(legendName), " = ", \
            histoList[histoNumber].Integral(0, histoList[histoNumber].FindBin(0.)-1)/histoList[histoNumber].Integral()
    
    if histoNumber == 0:
        histoList[histoNumber].Draw()
    else:
        histoList[histoNumber].Draw("same")

    if histoNumber == len(histoList)-1:
        legend.Draw("same")

pad2.cd()
gPad.SetLogy(0)
ratioList=[]
for histoNumber in range(len(histoList)):
    ratio = createRatio(histoList[histoNumber], histoList[histoNumber], "Top quark p_{T}", "#frac{Pythia 8}{Herwig 7}")
    ratio.SetDirectory(0)
    ratio = setHistoProperties(ratio, xTitle, yTitleRatio , ratioMin, ratioMax)

    ratioList.append(ratio)


for ratioNumber in range(len(ratioList)):
    #ratioList[ratioNumber].Print("all")
    # printUncertaintiesToLatex(ratioList[ratioNumber],ratioList[0],listOfLegendNames[ratioNumber])
    if ratioNumber == 0:
        ratioList[ratioNumber].Draw()
    else:
        ratioList[ratioNumber].Draw("same")


if not os.path.exists(outputPath):
  
  # Create a new directory because it does not exist 
  os.makedirs(outputPath)
  print("The new directory is created!")

canvas.Update()
string_normalize = ""
if normalize:
    string_normalize = "_normalized"

canvas.SaveAs(outputPath+"/"+variable+string_normalize+".pdf")

