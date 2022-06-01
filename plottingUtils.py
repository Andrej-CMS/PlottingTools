import sys
import ROOT
import ROOT.TH1D as TH1D
import ROOT.TH2D as TH2D
import ROOT.TGraph2D as TGraph2D
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TExec,TList,TString,TF1
from ROOT import kBlack, kBlue, kRed
from ROOT import TFile, gStyle,TGraphAsymmErrors, TLine, TGraph, TH1
import math
import array
import glob
import subprocess
import os
import re
import numpy as np
import pandas as pd



def createCanvasPads():
    c = TCanvas("c", "canvas", 1024, 1024)
    # Upper histogram plot is pad1
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    #pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.4)
    #pad2.SetGridx()
    pad2.Draw()
    return c, pad1, pad2

def createRatio(h1, h2, xTitle,yTitle):
    h3 = h1.Clone("ratio "+h1.GetTitle())
    h3.SetLineColor(h1.GetLineColor())
    h3.SetMarkerStyle(1)
    h3.SetTitle("")
    h3.SetMinimum(0.75)
    h3.SetMaximum(1.25)
    # Set up plot for markers and errors
    h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h2)
    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetTitle(yTitle)
    y.SetNdivisions(505)
    y.SetTitleSize(40)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.55)
    y.SetLabelFont(43)
    y.SetLabelSize(15)
    # Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetTitle(xTitle)
    #x.SetTitle("m[GeV]")
    #x.SetTitle("")
    x.SetTitleSize(40)
    x.SetTitleFont(43)
    x.SetTitleOffset(4.0)
    x.SetLabelFont(43)
    x.SetLabelSize(15)
    return h3

def getPadAxisHisto(pad):

    # Checking all objects existing in the pad
    padList = pad.GetListOfPrimitives()
    for objId in range(padList.GetEntries()):
        obj = padList.At(objId)
        if(obj.InheritsFrom("TH1")):
            print "object is TH1"
            return obj
        elif(obj.InheritsFrom("TGraph")): 
            print "object is TGraph"
            return obj.GetHistogram()


    print "Error in getPadAxisHisto: Couldn't get histogram for the axis"
    print "Histogram needs to be drawn on canvas/pad before making ratio pad"
    sys.exit()

def setHistoProperties(histo, xTitle, yTitle, minimum = 0.85, maximum = 1.15):
    h3 = histo.Clone("ratio "+histo.GetTitle())
    h3.SetTitle("")
    h3.SetMinimum(minimum)
    h3.SetMaximum(maximum)
    # Set up plot for markers and errors
    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetTitle(yTitle)
    y.CenterTitle()
    y.SetNdivisions(505)
    y.SetTitleSize(30)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.4)
    y.SetLabelFont(43)
    y.SetLabelSize(25)
    # Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetTitle(xTitle)
    x.SetTitleSize(40)
    x.SetTitleFont(43)
    x.SetTitleOffset(4.0)
    x.SetLabelFont(43)
    x.SetLabelSize(25)

    return h3

def getListOfVariableProperties(key):
    #Dictionary variable names corresponding to properties of the histogram
    # retuns a list of properties
    # List is structured as follows
    # [fixed order name, X-axis title, Y-axis title, nBins, [xBinsLow, xBinsHigh], [ratioMin, ratioMax], histo_totalXSname]
    variationDict = {
    "MC_j_afterFSR_pt" : ["pt1stjet_pt(ttb)>_0", "p_{T,1st jet}[GeV]","#frac{d#sigma }{dp_{T,1st jet} }", 
    						4, [20.,100.], [0., 2.01],
    						"Hardest additional parton jet", "total_pt(ttb)>_0", "truth"],

    "MC_sj_afterFSR_pt" : ["pt_jet_in_rho(ttbar+jet)>_0", "p_{T,#rho jet}[GeV]","#frac{d#sigma }{dp_{T,#rho jet} }", 
                            4, [20.,100.], [0., 2.01],
                            "#rho parton jet", "total_pt(ttb)>_0", "truth"],

    "MC_t_afterFSR_pt" : ["t_pt_pt(ttb)>_0", "p_{T,Top}[GeV]","#frac{d#sigma }{dp_{T,Top} }", 
    						13, [0.,1300.], [0., 2.01],
    						"Parton level top quark", "total_pt(ttb)>_0", "truth"],

    "MC_t_status62_pt" : ["t_pt_pt(ttb)>_0", "p_{T,Top}[GeV]","#frac{d#sigma }{dp_{T,Top} }", 
                            12, [20.,500.], [0., 2.01],
                            "Parton level top quark, Status 62", "total_pt(ttb)>_0", "truth"],

    "MC_rho_afterFSR" : ["rho_fb(ttbar+jet)>_0", "#rho","#frac{d#sigma }{d#rho }", 
                            4, [0.,1.], [0., 2.01],
                            "#rho distribution", "total_pt(ttb)>_0", "truth"],

    "MC_rho_sj_afterFSR" : ["rho_fb(ttbar+jet)>_0", "#rho","#frac{d#sigma }{d#rho }", 
                            4, [0.,1.], [0., 2.01],
                            "#rho distribution", "total_pt(ttb)>_0", "truth"],

    "jet_pt" : ["None", "p_{T,jet}[GeV]","#frac{d#sigma }{dp_{T,jet} }", 
                            15, [0.,1500.], [0., 2.01],
                            "Particle Level Jet p_{T}", "total_pt(ttb)>_0", "particleLevel"],
    "jet_HT" : ["None", "HT[GeV]","#frac{d#sigma }{dp_{HT} }", 
                            29, [250.,3150.], [0., 2.01],
                            "Particle Level HT", "total_pt(ttb)>_0", "particleLevel"],     
    "weight_mc" : ["None", "weight_{MC}","Events", 
                            120, [-2000.,10000.], [0., 1.99],
                            "MC weight", "total_pt(ttb)>_0", "truth"],                    
    }

    return variationDict.get(key)

    # End of getListOfVariableProperties

def getLegendNames(key):
    
    legendDict = {
        "run_HVQ_ATLAS_Standard"        : "Default",
        "run_HVQ_Bsup5_doubleFSR0"      : "Bsup 5, dFSR 0",
        "run_HVQ_Bsup5_doubleFSR1"      : "Bsup 5, dFSR 1",
        "run_HVQ_Bsup10_doubleFSR0"     : "Bsup 10, dFSR 0",
        "run_HVQ_Bsup10_doubleFSR1"     : "Bsup 10, dFSR 1",
        "run_HVQ_Bsup50_doubleFSR0"     : "Bsup 50, dFSR 0",
        "run_HVQ_Bsup50_doubleFSR1"     : "Bsup 50, dFSR 1",
        "run_HVQ_Bsup100_doubleFSR0"    : "Bsup 100, dFSR 0",
        "run_HVQ_Bsup100_doubleFSR1"    : "Bsup 100, dFSR 1",
        "run_HVQ_Bsup200_doubleFSR0"    : "Bsup 200, dFSR 0",
        "run_HVQ_Bsup200_doubleFSR1"    : "Bsup 200, dFSR 1",
        "run_HVQ_Bsup250_doubleFSR0"    : "Bsup 250, dFSR 0",
        "run_HVQ_Bsup250_doubleFSR1"    : "Bsup 250, dFSR 1",
        "run_HVQ_Bsup300_doubleFSR0"    : "Bsup 300, dFSR 0",
        "run_HVQ_Bsup300_doubleFSR1"    : "Bsup 300, dFSR 1",
        "run_HVQ_Bsup500_doubleFSR0"    : "Bsup 500, dFSR 0",
        "run_HVQ_Bsup500_doubleFSR1"    : "Bsup 500, dFSR 1",
    }

    return legendDict.get(key)


# function takes ratio histogram and prints uncertainties to a latex table
# given a second (nominal) ratio histogram, the function prints change in uncertainties
# compared to the nominal 

def printUncertaintiesToLatex(ratio, nominal, histoTitle):

    #print header

    print "\\begin{table}[htb]\\tiny"
    print "\\centering"
    print "\\caption{Change in statistical uncertainty $\\frac{\\sigma_{new}}{\\sigma_{def}}$ for ", histoTitle, " compared to the default ATLAS HVQ settings. }"
    #\resizebox{\textwidth}{!}{\begin{tabular}{lccccc}
    print "\t\\begin{tabular}{lc}"
    print "\t\t\\toprule"
    print "\t\tBin edges \t & \t $\\frac{\\sigma_{new}}{\\sigma_{def}}$ \\\\ \\midrule"

    # loop over all bins and calculate the difference in uncertainties
    # assumes same binning in ratio and nominal histo

    sumOfDifferences                = 0.
    sumOfDifferencesTo500           = 0.
    numBinsBelow500                 = 0.
    sumOfDifferences500To1000       = 0.
    numBins500To1000                = 0.
    sumOfDifferences1000andhigher   = 0.
    numBins1000andhigher            = 0.

    # lists for standard deviation
    sigmaTotal          = []
    sigmaBelow500       = []
    sigma500To1000      = []
    sigma1000andhigher  = []

    for xBin in range(ratio.GetNbinsX()):
        xBinLowerEdge = ratio.GetBinLowEdge(xBin+1)
        xBinUpperEdge = ratio.GetXaxis().GetBinUpEdge(xBin+1)

        uncRatio  = ratio.GetBinError(xBin+1)
        uncNominal= nominal.GetBinError(xBin+1)
        if histoTitle == "run_HVQ_ATLAS_Standard":
            differnceInUncertainties = uncRatio
        else:
            differnceInUncertainties = uncRatio/uncNominal
        if xBinLowerEdge<500.:
            numBinsBelow500 += 1.
            sumOfDifferencesTo500 += differnceInUncertainties
            sigmaBelow500.append(differnceInUncertainties)
        elif xBinLowerEdge < 1000.:
            numBins500To1000 += 1.
            sumOfDifferences500To1000 += differnceInUncertainties
            sigma500To1000.append(differnceInUncertainties)
        else:
            numBins1000andhigher += 1.
            sumOfDifferences1000andhigher+= differnceInUncertainties
            sigma1000andhigher.append(differnceInUncertainties)
        sumOfDifferences += differnceInUncertainties
        sigmaTotal.append(differnceInUncertainties)

        print "\t\t$[",xBinLowerEdge, ",", xBinUpperEdge ,"]$", "\t & \t" , differnceInUncertainties , "\\\\"

    # print end of table
    print "\t\t\\bottomrule"
    print "\t\\end{tabular}"
    print "\\end{table} \n"

    meanTotal = sumOfDifferences/ratio.GetNbinsX()
    meanBelow500 = sumOfDifferencesTo500/numBinsBelow500
    mean500To1000      = sumOfDifferences500To1000/numBins500To1000
    mean1000andhigher  = sumOfDifferences1000andhigher/numBins1000andhigher

    # Calculate Standard Deviation

    stdDevTotal         = 0.
    stdDevBelow500      = 0.
    stdDev500To1000     = 0.
    stdDev1000andhigher = 0.
    for xBin in range(ratio.GetNbinsX()):
        stdDevTotal         += (sigmaTotal[xBin] - meanTotal)**2

    for xBin in range(len(sigmaBelow500)):
        stdDevBelow500      += (sigmaBelow500[xBin] - meanBelow500)**2

    for xBin in range(len(sigma500To1000)):
        stdDev500To1000     += (sigma500To1000[xBin] - mean500To1000)**2

    for xBin in range(len(sigma1000andhigher)):
        stdDev1000andhigher += (sigma1000andhigher[xBin] - mean1000andhigher)**2

    print "Average change in uncertainty  $= ", meanTotal, " \\pm ", math.sqrt(stdDevTotal/ratio.GetNbinsX()), "$"
    print "Average below 500   $= ", meanBelow500        , " \\pm ", math.sqrt(stdDevBelow500/numBinsBelow500), "$"
    print "Average 500-1000    $= ", mean500To1000       , " \\pm ", math.sqrt(stdDev500To1000/numBins500To1000), "$"
    print "Average above 1000  $= ", mean1000andhigher   , " \\pm ", math.sqrt(stdDev1000andhigher/numBins1000andhigher), "$"
    print "-------------------------------------------------"
# END OF printUncertaintiesToLatex
