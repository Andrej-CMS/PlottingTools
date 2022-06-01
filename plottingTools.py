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
#from scipy.special import zeta, polygamma, factorial

#sys.path.insert(0, './rundec-python/build/lib.linux-x86_64-2.7')
#import rundec


def alphaS(alphaSatZMass,scale):
    nFlav=5
    xeta3=zeta(3)
    qOoverLambda=scale/0.0355
    beta0=11-2/3.*nFlav
    beta1=102-38/3*nFlav
    beta2=2857/2.-5033/18.*nFlav+325/54.*nFlav**2
    beta3=(149753/6.+3564*xeta3)-(1078361./162+6508/27.*xeta3)*nFlav+(50065.0/162+6472/81.0*xeta3)*nFlav**2+1093.0/729*nFlav**3

    alpha_S=4*math.pi/(beta0*math.log(qOoverLambda**2))*(1-(beta1*math.log(math.log(qOoverLambda**2)))/(beta0**2*math.log(qOoverLambda**2)) \
        +beta1**2/(beta0**4*math.log(qOoverLambda**2)**2)*((math.log(math.log(qOoverLambda**2))**2)-math.log(math.log(qOoverLambda**2))-1-beta2*beta1/(beta1**2))  \
        +beta1**3/(beta0**6*math.log(qOoverLambda**2)**3)*(5/2.0*math.log(math.log(qOoverLambda**2))**2-math.log(qOoverLambda**2)**3   \
        +2*math.log(math.log(qOoverLambda**2))-1/2.0-3*beta2*beta0/(beta1**2)*math.log(math.log(qOoverLambda**2))+beta3*beta0**2/(2*beta1**3)))
    return alpha_S


def poleMass(runningMass, scale):
    alpha_S=alphaS(0.118,scale)
    nFlav=5 #number of active flavors
    logarithm=math.log(float(scale**2)/runningMass**2)
    xeta3=zeta(3)
    d1=4/3.+logarithm
    d2=307/32.0+(math.pi**2)/3+((math.pi)**2)/9*math.log(2.)-1/6.*xeta3+509/72.*logarithm+47/24.*(logarithm**2)-nFlav*(71/144.+(math.pi**2)/18.+13/36.*logarithm+1/12.*logarithm**2)
    #now combine for top pole mass calculation
    return runningMass*(1+alpha_S/math.pi*d1+(alpha_S/math.pi)**2*d2)

def getRunningMass(targetPoleMass,scale_factor,runningMass_v=[150.,180.]):
    precision=0.004
    HiggsMass=125.0
 
    lowerValue=poleMass(runningMass_v[0], scale_factor*(2*runningMass_v[0]+HiggsMass))
    upperValue=poleMass(runningMass_v[1], scale_factor*(2*runningMass_v[1]+HiggsMass))

    if abs(targetPoleMass-lowerValue)<precision:
        return runningMass_v[0],scale_factor*(2*runningMass_v[0]+HiggsMass)
    elif abs(upperValue-targetPoleMass)<precision:
        #print float(runningMass_v[1]) , scale_factor*(2*runningMass_v[1]+HiggsMass)
        return float(runningMass_v[1]) , scale_factor*(2*runningMass_v[1]+HiggsMass)

    elif lowerValue<targetPoleMass and upperValue>targetPoleMass:
        lowerBound=runningMass_v[0]
        upperBound=(runningMass_v[1]+runningMass_v[0])/2

    elif lowerValue>targetPoleMass :
        lowerBound=runningMass_v[0]-(runningMass_v[1]-runningMass_v[0])
        upperBound=runningMass_v[0]
    elif upperValue<targetPoleMass:
        lowerBound=runningMass_v[1]
        upperBound=runningMass_v[1]+(runningMass_v[1]-runningMass_v[0])
        
    return getRunningMass(targetPoleMass,scale_factor,[lowerBound,upperBound])
    



def createCanvasPads():
    c = TCanvas("c", "canvas", 800, 800)
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

def createCanvas3Pads():
    c = TCanvas("c", "canvas", 1000, 1000)

    # Upper histogram plot is pad1
    pad1 = TPad("pad1", "pad1", 0, 0.6, 1, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot

    #pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2

    c.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.4, 1, 0.6)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0)
    
    #pad2.SetGridx()
    pad2.Draw()

    c.cd()
    pad3 = TPad("pad3", "pad3", 0, 0.1, 1, 0.4)
    pad3.SetTopMargin(0)  # joins upper and lower plot
    pad3.SetBottomMargin(0.4)
    
    #pad2.SetGridx()
    pad3.Draw()
    return c, pad1, pad2, pad3


def createRatio(h1, h2, xTitle,yTitle):
    h3 = h1.Clone("h3")
    h3.SetLineColor(ROOT.kRed)
    h3.SetMarkerStyle(1)
    h3.SetTitle("")
    h3.SetMinimum(0.85)
    h3.SetMaximum(1.15)
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

def setHistoProperties(histo, xTitle, yTitle, minimum = 0.85, maximum = 1.15):
    h3 = histo.Clone("h3")
    h3.SetTitle("")
    h3.SetMinimum(minimum)
    h3.SetMaximum(maximum)
    # Set up plot for markers and errors
    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetTitle(yTitle)
    y.CenterTitle()
    y.SetNdivisions(505)
    y.SetTitleSize(40)
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

def makeHistogramFromMG5File(fileName, variableName, variationNumber=1):
# simple code to turn MG5 output files into arrays of floats that can be fed into a histogram for plotting
    print "Variable in makeHistogramFromMG5File " , variableName
    found = False
    skipLines = 0
    variation=1 #variable needed to skip histograms in top drawer file until the program gets muR=muF variation
    collumns = []
    binCenters = []
    binEntries = []
    binUncertainties = []
    with open(fileName) as openFile:
        for line in openFile:
            #find variable
            if variableName in line:
                found=True
                if variation < 2*variationNumber-1: #top drawer has two entried with the name per histogram
                    found=False
                    variation+=1
                
                
            if not found:
                continue
            #skip not required lines in MG5 file
            if skipLines<7:
                skipLines+=1
                continue
            #end of histogram
            if "HIST SOLID" in line:
                break
            collumns.append(line.split())
    #change numbers into scientific format
    collumns = [[x.replace('D','E') for x in l] for l in collumns]
    for i in range(len(collumns)):
        binCenters.append(float(collumns[i][0]))
        binEntries.append(float(collumns[i][1]))
        binUncertainties.append(float(collumns[i][2]))
    return binCenters, binEntries, binUncertainties

def fillAndReturnHistogram(histoName,binCenters,binEntries,binUncertainties):
    #TH1D needs array
    binLowEdges = makeBinLowEdges(binCenters)
    histo = TH1D(histoName, histoName, len(binCenters), binLowEdges)
    #fill histo
    for iterator in range(histo.GetNbinsX()):
        histo.AddBinContent(histo.FindBin(binCenters[iterator]),binEntries[iterator])
        histo.SetBinError(histo.FindBin(binCenters[iterator]),binUncertainties[iterator])  
    #Checked: histograms are the same (binCenter,Value,Error) as in .top files from MG5
    return histo
def getBinLowEdges(histo):
    xAxislowEdges= array.array('d')
    for bin in range(histo.GetNbinsX()+1):
        xAxislowEdges.append(histo.GetXaxis().GetBinLowEdge(bin+1))
    
    #checked that xAxislowEdges and  histogram edges are the same
    return xAxislowEdges

def make1DHistoList(folderName, fileName, variableInMadgraph,massRange,scale=-1.,variation=-1.):
    histoList1D= []
#    path = folderName+fileName
#    for filepath in glob.glob(path):
#        mass = float(filepath[filepath.find("_M_")+len("_M_"):filepath.rfind("_Mu_")])
#        print filepath, mass
#        binCenters,binEntries,binUncertainties = makeHistogramFromMG5File(filepath, variableInMadgraph)
#        histoList1D.append(fillAndReturnHistogram(str(mass),binCenters,binEntries,binUncertainties))
    if variation <0:
        
        if scale < 0.:
            for mass in massRange:
                #check if MG5 output exists
                try:
                    open(folderName+str(mass)+fileName)
                except:
                    print "WARNING: ", folderName+str(mass)+fileName , " doesn't exist!!! "
                    continue
                binCenters,binEntries,binUncertainties = makeHistogramFromMG5File(folderName+str(mass)+fileName, variableInMadgraph)
                histoList1D.append(fillAndReturnHistogram(str(mass),binCenters,binEntries,binUncertainties))
        else:

            for mass in massRange:
                #check if MG5 output exists
                filePath=folderName+str('{0:g}'.format(mass))+"_Mu_"+str(scale).rstrip("0")+fileName
                try:
                    print "Open File " , filePath
                    open(filePath)
                except:
                    print "WARNING: ", filePath , " doesn't exist!!! "
                    continue
                binCenters,binEntries,binUncertainties = makeHistogramFromMG5File(filePath, variableInMadgraph)
                histoList1D.append(fillAndReturnHistogram(str(mass),binCenters,binEntries,binUncertainties))
    else:
        # get MuR, muF variation from Top Drawer File
        variationNumber = GetVariationEntryInFile("{0:.3f}".format(variation))
        for mass in massRange:
            filePath=folderName+str('{0:g}'.format(mass))+"_Mu_"+str(scale).rstrip("0")+fileName
            variableInMadgraph=variableInMadgraph+" dyn=   0 muR= "+"{0:.3f}".format(variation)#+" muF= "+"{0:.3f}".format(variation)
            #print "variableInMadgraph" , variableInMadgraph
            
            try:
                print "Open File " , filePath
                open(filePath)
            except:
                print "WARNING: ", filePath , " doesn't exist!!! "
                continue
            binCenters,binEntries,binUncertainties = makeHistogramFromMG5File(filePath, variableInMadgraph,variationNumber)
            histoList1D.append(fillAndReturnHistogram(str(mass)+" variation "+str(variation),binCenters,binEntries,binUncertainties))
    #Histograms filled correctly from top file. checked!
    #histoList1D[0].Print("all")

    return histoList1D
    

def make2DHistfromListOf1DHisto(histoList1D = [] , massRange = [] ,name='',title=''):
    massLowEdges = makeBinLowEdges(massRange)
    yAxisLowEdges = getBinLowEdges(histoList1D[0])
    histoTest= TH2D (name, title, len(massRange), massLowEdges, histoList1D[0].GetNbinsX(), yAxisLowEdges)

    histo2D = TGraph2D()
    #histoTest = TH2D(name, title, len(massRange)-1, massRange, len(histoList1D)-1, histoList1D)
    pointIt=0
    for massIt in range(len(massRange)):
        currentHisto=histoList1D[massIt]
        print currentHisto.GetName()
        for binIt in range(histoList1D[massIt].GetNbinsX()+2):
            histo2D.SetPoint(pointIt,massRange[massIt],currentHisto.GetXaxis().GetBinCenter(binIt+1),currentHisto.GetBinContent(binIt+1))
            pointIt+=1
            histoTest.Fill (massRange[massIt], currentHisto.GetXaxis().GetBinCenter(binIt+1), currentHisto.GetBinContent(binIt+1))
            histoTest.SetBinError(massIt+1,binIt+1,currentHisto.GetBinError(binIt+1))
        #cehcked that 2D histo has the correct entries for [mass,HistoBin] tuple
    return histo2D,histoTest

def getBornXSDerivativeFrom2DHisto(pathToFile, distribution, mass, bin_center, xBins=[], yBins=[]):
#    print "getBornXSDerivativeFrom2DHisto, File = ", pathToFile , " Distribution " ,distribution
    inputFile = TFile(pathToFile)
    histo = TH2D(inputFile.Get(distribution))
    #histo.Sumw2()

    if len(xBins) and len(yBins):
        histo = rebin2DVariableBinning(histo, xBins, yBins)
    
    #expects X-axis:Mass, Y-axis:distribution that one wants to look at (e.g. pT of the Higgs boson)
    massBin = histo.GetXaxis().FindBin(mass)
    valueBin = histo.GetYaxis().FindBin(bin_center)
    bornDerivate=0.
    if massBin == histo.GetNbinsX():
        bornDerivate = (histo.GetBinContent(massBin-1, valueBin)-histo.GetBinContent(massBin, valueBin))/(histo.GetXaxis().GetBinCenter(massBin-1)-histo.GetXaxis().GetBinCenter(massBin))
    else:
        bornDerivate = (histo.GetBinContent(massBin+1, valueBin)-histo.GetBinContent(massBin, valueBin))/(histo.GetXaxis().GetBinCenter(massBin+1)-histo.GetXaxis().GetBinCenter(massBin))

    return bornDerivate    


def makeListOfMassDependenceHistos(histoList1D = [] , massRange = []):
    massDependenceHistos = []
    binLowEdges = makeBinLowEdges(massRange)
    #Loop over all Masspoints and Histograms to create histograms that show 
    #the mass dependence of the cross section in each bin of the differential distribution (pT, eta, etc.)
    for theBin in range(histoList1D[0].GetNbinsX()+2):
        massHisto = TH1D (str(histoList1D[0].GetBinCenter(theBin)), str(histoList1D[0].GetBinCenter(theBin)), len(massRange), binLowEdges)
        for histoID, histo in  enumerate(histoList1D):
            massHisto.SetBinContent(histoID,histo.GetBinContent(theBin))
            massHisto.SetBinError(histoID,histo.GetBinError(theBin))
          #  print histoID, histo.GetBinContent(theBin), theBin,histo.GetBinError(theBin)
          ## set plot style
            y = massHisto.GetYaxis()
            y.SetTitle("#frac{d#sigma}{dM_{Top}} ")
            y.SetNdivisions(505)
            y.SetTitleSize(20)
            y.SetTitleFont(43)
            y.SetTitleOffset(1.55)
            y.SetLabelFont(43)
            y.SetLabelSize(15)
            # Adjust x-axis settings
            x = massHisto.GetXaxis()
            x.SetTitle("p_{t}[GeV]")
            #x.SetTitle("m[GeV]")
            #x.SetTitle("")
            x.SetTitleSize(20)
            x.SetTitleFont(43)
            x.SetTitleOffset(4.0)
            x.SetLabelFont(43)
            x.SetLabelSize(15)
            massHisto.SetMinimum(0.001)
         
        massDependenceHistos.append(massHisto)

    return massDependenceHistos

def makeBinLowEdges(bins = []):
    #tested, the low bins are correct!
    #assumes the bins to be sorted
    binLowEdges = array.array('d')
    arraysize=len(bins)
    for binN in range(arraysize):
        if binN == 0:
            binLowEdges.append(bins[binN]-(float(bins[binN+1])-float(bins[binN]))/2)
        else:
            #binLowEdges.append(bins[binN]-(float(bins[binN])-float(bins[binN-1]))/2)
            binLowEdges.append((float(bins[binN])+float(bins[binN-1]))/2)
        #add upper bin edge for last bin
        if binN == arraysize-1:
            binLowEdges.append((float(bins[binN])-float(bins[binN-1]))/2+float(bins[binN]))
    return binLowEdges

def fitListofHistos(listOfHistos):
    for histo in listOfHistos:
        histo.Fit("pol1","F")
        #histo.SetStats(0)
        gStyle.SetOptFit(1111)
        # Set stat options
        gStyle.SetStatY(0.5)                
        # Set y-position (fraction of pad size)
        gStyle.SetStatX(0.7)                
        # Set x-position (fraction of pad size)
        #gStyle.SetStatW(0.4)                
        # Set width of stat-box (fraction of pad size)
        #gStyle.SetStatH(0.2)                
        # Set height of stat-box (fraction of pad size)

    return listOfHistos

def calculateGradient(histo2D):
    #histo2D=histo2DTGraph.GetHistogram()
    #histo2D.Draw("surf1") 
    debug = False
    
    xArray = np.zeros(histo2D.GetXaxis().GetNbins())
    yArray = np.zeros(histo2D.GetYaxis().GetNbins())
    if debug:
        print len(xArray),len(yArray)
    for xBin in range(histo2D.GetXaxis().GetNbins()):
        xArray[xBin]=histo2D.GetXaxis().GetBinCenter(xBin)
    for yBin in range(histo2D.GetYaxis().GetNbins()):
        yArray[yBin]=histo2D.GetYaxis().GetBinCenter(yBin)
        
    X, Y = np.meshgrid(xArray, yArray)
    zs = np.array([histo2D.GetBinContent(histo2D.GetXaxis().FindBin(xArray),histo2D.GetXaxis().FindBin(yArray)) for xArray,yArray in zip(np.ravel(X), np.ravel(Y))])
    #print zs
    Z = zs.reshape(X.shape)

    gx,gy = np.gradient(Z,0.5,10)
    return gx,gy


def calcuateDiffXSRunningMass(nloXSHisto, pdf, runningMass, scale, variableInMadgraph,xBins=[],yBins=[]):
    # Get derrivative of born XS
    ## read from root file containing (2D) histogram with mass dependence
    ##
    debug = False
    crd = rundec.CRunDec()
    runningMassHisto = nloXSHisto[0].Clone("runningMassHisto")
    scaleFactor = float(scale)/(2*runningMass+125.0)
    list1=[0.125, 0.167, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0]
    scaleFactor=min(list1, key=lambda x:abs(x-scaleFactor)) # choose closest value to the one in list (avoid rounding errors)
    print "scale in calcuateDiffXSRunningMass " , scale , " scale factor ", scaleFactor
    #alphas=alphaS(0.118,scale)
    alphas=crd.AlphasExact(0.118, 91.1876, scale, 5, 2)
    diffBornXSList=[]
    for currentBin in range(nloXSHisto[0].GetNbinsX()+2):

        diffBornXS= getBornXSDerivativeFrom2DHisto("MassdependenceRootFiles_PreciseGrids/ttH_BORN_MassDependence_"+pdf+"_Mu_"+str(scaleFactor)+".root", variableInMadgraph, runningMass, nloXSHisto[0].GetBinCenter(currentBin),xBins,yBins)
        diffBornXSList.append(diffBornXS)
        logarithm=math.log(float(scale)**2/runningMass**2)
        d1=4/3.+logarithm
        if debug:
            print "########## Values for XS Combination ##########"
            print "Born derivative, correction , XS value NLO in bin  : ", diffBornXS, alphas/math.pi*d1*runningMass*diffBornXS , nloXSHisto[0].GetBinContent(currentBin)
        newValue = nloXSHisto[0].GetBinContent(currentBin)+(alphas/math.pi)*d1*runningMass*diffBornXS
        runningMassHisto.SetBinContent(currentBin, newValue)
        runningMassHisto.SetBinError(currentBin,nloXSHisto[0].GetBinError(currentBin))
    if debug:
        print diffBornXSList
    return runningMassHisto

def makeDiffXSPlotUncertainties(inputHistos):
    nBins = inputHistos[0].GetNbinsX()
    x = array.array('d')
    exl = array.array('d')
    exh = array.array('d')
    y = array.array('d')
    eyl = array.array('d')
    eyh = array.array('d')

    for currentBin in range(nBins+1):
        x.append(inputHistos[0].GetBinCenter(currentBin))
        y.append(inputHistos[0].GetBinContent(currentBin))

        #no error in X 
        exl.append(0.0)
        exh.append(0.0)

        if inputHistos[1].GetBinContent(currentBin)>inputHistos[0].GetBinContent(currentBin):
            eyh.append(abs(inputHistos[1].GetBinContent(currentBin)-inputHistos[0].GetBinContent(currentBin)))
            eyl.append(abs(inputHistos[0].GetBinContent(currentBin)-inputHistos[2].GetBinContent(currentBin)))
        else:
            eyh.append(abs(inputHistos[2].GetBinContent(currentBin)-inputHistos[0].GetBinContent(currentBin)))
            eyl.append(abs(inputHistos[0].GetBinContent(currentBin)-inputHistos[1].GetBinContent(currentBin)))
    
    assymErrorGraph_XS = TGraphAsymmErrors(nBins,x,y,exl,exh,eyl,eyh)

    return assymErrorGraph_XS

def makeDiffXSGraphError(inputHistos,xError=False,shiftValues=False):
    debug = False
    #list [[up,down]] with relative uncertainties (normalized to the nominal value)
    errors, errorsUp, errorsDown = extractUncertaintiesFromHistos(inputHistos)

    if debug:
        print "Uncertainties in makeDiffXSGraphError "
        print errors
        print "####Errors UP##########"
        print errorsUp
    nBins = len(errors)

    if nBins != inputHistos[0].GetNbinsX():
        print "Warning: Number of bins in errors and histograms are different"
        print "Setting all other errors to 0"

    graph = TGraphAsymmErrors(inputHistos[0])
    binValues = array.array('d',[])
    binCenters = array.array('d',[])
    exl = array.array('d',[])
    exh = array.array('d',[])
    #set uncertainties
    for xBin in range(nBins):
        binValue = inputHistos[0].GetBinContent(xBin+1)
        binValues.append(binValue)
        binCenters.append(inputHistos[0].GetBinCenter(xBin+1))
        exl.append(0.)
        exh.append(0.)

        
        #needed for x-axis range of each bin
        binCenter = inputHistos[0].GetBinCenter(xBin+1)
        #upper bin edge
        binXUp = inputHistos[0].GetBinLowEdge(xBin+1)
        #lower bin edge
        binXDown = inputHistos[0].GetBinLowEdge(xBin)
        if shiftValues:
            graph.SetPoint(xBin, binCenter+0.2*inputHistos[0].GetBinWidth(xBin+1), inputHistos[0].GetBinContent(xBin+1))


        if(np.sign(errors[xBin][0]) != np.sign(errors[xBin][1])):
            graph.SetPointEYhigh(xBin, abs(errors[xBin][0])*binValue)
            graph.SetPointEYlow(xBin, abs(errors[xBin][1])*binValue)

        elif(np.sign(errors[xBin][0])>0):
            if abs(errors[xBin][0])>= abs(errors[xBin][1]):
                graph.SetPointEYhigh(xBin, abs(errors[xBin][0])*binValue)
                graph.SetPointEYlow(xBin, 0.)
            else:
                graph.SetPointEYhigh(xBin, abs(errors[xBin][1])*binValue)
                graph.SetPointEYlow(xBin, 0.)

        elif(np.sign(errors[xBin][0])<0):
            if abs(errors[xBin][0])>= abs(errors[xBin][1]):
                graph.SetPointEYhigh(xBin, 0.)
                graph.SetPointEYlow(xBin, abs(errors[xBin][0])*binValue)
            else:
                graph.SetPointEYhigh(xBin, 0.)
                graph.SetPointEYlow(xBin, abs(errors[xBin][1])*binValue)

        if(xError == False):
            graph.SetPointEXhigh(xBin, 0.);
            graph.SetPointEXlow(xBin, 0.);
        else:
            graph.SetPointEXhigh(xBin, inputHistos[0].GetBinWidth(xBin)/2);
            graph.SetPointEXlow(xBin, inputHistos[0].GetBinWidth(xBin)/2);

        # Latex print central values and uncertainties for publication


#    graph = TGraphAsymmErrors(nBins,binCenters,binValues,exl,exh,errorsDown,errorsUp)

    return graph

def extractUncertaintiesFromHistos(inputHistos, symmetricErrors = False):
        nBins = inputHistos[0].GetNbinsX()
        #relative variations in each bin [[up,down]]
        errors = [[0.,0.]]*nBins
        errorsUp = array.array('d',[])
        errorsDown = array.array('d',[])

        for iBin in range(nBins):
            #leave out underflow bin
            nominalValue = inputHistos[0].GetBinContent(iBin+1)
            histUpValue = inputHistos[1].GetBinContent(iBin+1)
            histDownValue = inputHistos[2].GetBinContent(iBin+1)
            relativeVarUp = 0.
            relativeVarDown = 0.
            if 0 == nominalValue:
                print "WARNING : nominalValue for histogram is 0"
                print "Setting errors to up/ down value"
                relativeVarUp = histUpValue
                relativeVarDown = histDownValue

            else:

                if (histUpValue >= nominalValue) and (histDownValue <= nominalValue):
                    errorsUp.append(histUpValue-nominalValue)
                    errorsDown.append(nominalValue-histDownValue)
                    relativeVarUp = (histUpValue-nominalValue)/nominalValue
                    relativeVarDown = (histDownValue-nominalValue)/nominalValue

                elif (histUpValue <= nominalValue) and (histDownValue >= nominalValue):
                    errorsUp.append(histDownValue-nominalValue)
                    errorsDown.append(nominalValue-histUpValue)
                    relativeVarUp = (histDownValue-nominalValue)/nominalValue
                    relativeVarDown = (histUpValue-nominalValue)/nominalValue

                elif (np.sign(histUpValue-nominalValue) == np.sign(histDownValue-nominalValue)):
                    #print "Warning in extractUncertaintiesFromHistos : Up and Down variations have same sign"
                    if histUpValue-nominalValue > 0:
                        relativeVarDown = 0.
                        errorsUp.append(abs(histUpValue-nominalValue))
                        errorsDown.append(0.)
                        if abs(histUpValue-nominalValue) > abs(histDownValue-nominalValue):
                            relativeVarUp = (histUpValue-nominalValue)/nominalValue
                        else:
                            relativeVarUp = (histDownValue-nominalValue)/nominalValue
                    else:
                        relativeVarUp = 0.
                        errorsUp.append(0.)
                        errorsDown.append(abs(histUpValue-nominalValue))
                        if abs(histUpValue-nominalValue) > abs(histDownValue-nominalValue):
                            relativeVarDown = (histUpValue-nominalValue)/nominalValue
                        else:
                            relativeVarDown = (histDownValue-nominalValue)/nominalValue

                        
                    
                    #relativeVarUp = (histUpValue-nominalValue)/nominalValue
                    #relativeVarDown = (histDownValue-nominalValue)/nominalValue

            errors[iBin] = [relativeVarUp,relativeVarDown]


        return errors , errorsUp, errorsDown

def setGraphStyle(graph, line, lineColor, lineWidth, 
                        marker, markerColor, markerSize,
                        fill, fillColor):

    if(line != -1):
        graph.SetLineStyle(line)
    if(lineColor != -1):
        graph.SetLineColor(lineColor)
    if(lineWidth != -1):
        graph.SetLineWidth(lineWidth)
    
    if(fill != -1):
        graph.SetFillStyle(fill)
    if(fillColor != -1):
        graph.SetFillColor(fillColor)
        graph.SetFillColorAlpha(fillColor, 0.35)
    
    if(marker != -1):
        graph.SetMarkerStyle(marker)
    if(markerColor != -1):
        graph.SetMarkerColor(markerColor)
    if(markerSize != -1):
        graph.SetMarkerSize(markerSize)


def makeRatioGraph(graph_nominator,graph_denominator,errorType,shiftValues=False):
    ratioGraph = TGraphAsymmErrors(graph_nominator.Clone());

    for iBin in range(ratioGraph.GetN()):
        nominator_v = ROOT.double(0)
        nominatorNext_v = ROOT.double(0)
        nominatorPrevious_v = ROOT.double(0)

        denominator_v = ROOT.double(0) 

        xValue = ROOT.double(0)
        xValueNextPoint = ROOT.double(0)
        xValuePreviousPoint = ROOT.double(0)

        ratio_v = ROOT.double(0)

        graph_nominator.GetPoint(iBin,xValue,nominator_v)
        graph_denominator.GetPoint(iBin,xValue,denominator_v)

        #caclulate bin width
        if iBin < ratioGraph.GetN()-1:
            graph_nominator.GetPoint(iBin+1,xValueNextPoint,nominatorNext_v)
            binWidth = abs(xValueNextPoint-xValue)/2.0
        else:
            graph_nominator.GetPoint(iBin-1,xValuePreviousPoint,nominatorPrevious_v)
            binWidth = abs(xValue-xValuePreviousPoint)/2.0

        if denominator_v > 0. :
            ratio_v = nominator_v / denominator_v
        else:
            ratio_v = 1.

        nominator_errorUP = graph_nominator.GetErrorYhigh(iBin)
        nominator_errorDOWN = graph_nominator.GetErrorYlow(iBin)

        denominator_errorUP = graph_denominator.GetErrorYhigh(iBin)
        denominator_errorDOWN = graph_denominator.GetErrorYlow(iBin)

        if denominator_v != 0.:
            if errorType == 1:
                ratio_errorUP = nominator_errorUP/denominator_v
                ratio_errorDOWN = nominator_errorDOWN/denominator_v
            else:
                ratio_errorUP = nominator_errorUP
                ratio_errorDOWN = nominator_errorDOWN
        else:
            ratio_errorUP = nominator_errorUP
            ratio_errorDOWN = nominator_errorDOWN

        if shiftValues:
            ratioGraph.SetPoint(iBin, xValue+0.2*binWidth, ratio_v)
            ratioGraph.SetPointEYhigh(iBin, ratio_errorUP)
            ratioGraph.SetPointEYlow(iBin, ratio_errorDOWN)
        else:
            ratioGraph.SetPoint(iBin, xValue, ratio_v)
            ratioGraph.SetPointEYhigh(iBin, ratio_errorUP)
            ratioGraph.SetPointEYlow(iBin, ratio_errorDOWN)
    
    return ratioGraph

def drawRatioPad(pad, yMin, yMax, title, fraction=0.36):
    #need first to get the histogram from the pad/canvas for axis to draw
    axisHisto = getPadAxisHisto(pad)

    # y:x size ratio for canvas
    canvAsym = (pad.GetY2() - pad.GetY1())/(pad.GetX2()-pad.GetX1())
    left  = pad.GetLeftMargin()
    right = pad.GetRightMargin()
    # change old pad
    pad.SetBottomMargin(fraction)
    pad.SetRightMargin(right)
    pad.SetLeftMargin(left)
    pad.SetBorderMode(0)
    pad.SetBorderSize(0)
    pad.SetFillColor(10)
    # create new pad for ratio plot
    rPad = TPad("rPad","",0,0,1,fraction+0.001)

    rPad.SetFillStyle(0)
    rPad.SetFillColor(0)

    rPad.SetBorderSize(0)
    rPad.SetBorderMode(0)
    rPad.Draw()
    rPad.cd()
    rPad.SetLogy(0)
#    rPad.SetTicky(1)
    # configure ratio plot
    scaleFactor = 1./(canvAsym*fraction)
    h_axis = axisHisto.Clone("h_axis")
    h_axis.SetStats(ROOT.kFALSE)
    h_axis.SetTitle("")
    h_axis.SetName("h_axis")
    h_axis.SetMaximum(yMax)
    h_axis.SetMinimum(0.8)
    # configure axis of the pad
    h_axis.GetXaxis().SetTitleSize(axisHisto.GetXaxis().GetTitleSize()*scaleFactor*1.3)
    h_axis.GetXaxis().SetTitleOffset(axisHisto.GetXaxis().GetTitleOffset()*0.9)
    h_axis.GetXaxis().SetLabelSize(axisHisto.GetXaxis().GetLabelSize()*scaleFactor*1.05)
    h_axis.GetXaxis().SetTitle(axisHisto.GetXaxis().GetTitle())
    h_axis.GetXaxis().SetNdivisions(axisHisto.GetNdivisions())
    h_axis.GetYaxis().CenterTitle()
    h_axis.GetYaxis().SetTitle(title)
    h_axis.GetYaxis().SetTitleSize(axisHisto.GetYaxis().GetTitleSize()*scaleFactor)
    h_axis.GetYaxis().SetTitleOffset(axisHisto.GetYaxis().GetTitleOffset()/scaleFactor)
    h_axis.GetYaxis().SetLabelSize(axisHisto.GetYaxis().GetLabelSize()*scaleFactor)
    h_axis.GetYaxis().SetLabelOffset(axisHisto.GetYaxis().GetLabelOffset()*3.3)
    h_axis.GetYaxis().SetTickLength(0.03)
    h_axis.GetYaxis().SetNdivisions(405)
    h_axis.GetXaxis().SetRange(axisHisto.GetXaxis().GetFirst(), axisHisto.GetXaxis().GetLast())
    # delete axis of initial plot
    axisHisto.GetXaxis().SetLabelSize(0)
    axisHisto.GetXaxis().SetTitleSize(0)
    #This is frustrating and stupid but apparently necessary...
    setex1 = TExec("setex1","gStyle.SetErrorX(0.5)") 
    setex1.Draw()
    setex2 = TExec("setex2","gStyle.SetErrorX(0.)") 
    setex2.Draw()
    h_axis.Draw("axis")
    rPad.SetTopMargin(0.0)
    rPad.SetBottomMargin(0.15*scaleFactor)
    rPad.SetRightMargin(right)
    pad.SetLeftMargin(left)
    pad.RedrawAxis()
    # draw grid
    rPad.SetGrid(0,1)
    rPad.cd()
    
    # draw a horizontal line on the pad
    xmin = h_axis.GetXaxis().GetXmin()
    xmax = h_axis.GetXaxis().GetXmax()
    height = TString("") 
    height += 1
    f = TF1("f", height, xmin, xmax)
    f.SetLineStyle(1)
    f.SetLineWidth(1)
    f.SetLineColor(kBlack)
    f.Draw("L same")
    
    
    return rPad

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


    
def writeHistosToFile(listOfHistos,outputFile):
    output= TFile(outputFile,"RECREATE")
    output.Write()
    output.Close()

def rebin2DVariableBinning(originalHisto, xBins, yBins):
    debug = False

    rebinnedHisto = TH2D ("Rebinned Histo", "Rebinned Histo", len(xBins)-1, xBins, len(yBins)-1, yBins)
    
    for xBin in range(1,rebinnedHisto.GetNbinsX()+2):
        #Get bin-edges X-axis
        xBinLowEdge=rebinnedHisto.GetXaxis().GetBinLowEdge(xBin)
        xBinUpperEdge=rebinnedHisto.GetXaxis().GetBinLowEdge(xBin)+rebinnedHisto.GetXaxis().GetBinWidth(xBin)
        
        #Find corresponding bins in original histogram
        orig_xBinLow = originalHisto.GetXaxis().FindBin(xBinLowEdge+0.001)
        orig_xBinUp = originalHisto.GetXaxis().FindBin(xBinUpperEdge-0.001)

        orig_xBinLowEdge = originalHisto.GetXaxis().GetBinLowEdge(orig_xBinLow)
        orig_xBinUpEdge = originalHisto.GetXaxis().GetBinLowEdge(orig_xBinUp)+originalHisto.GetXaxis().GetBinWidth(orig_xBinUp)
        #is correct
        if debug:
            print('New Bins [{:f}-{:f}] \t Original Bins [{:f}-{:f}]'.format(xBinLowEdge, xBinUpperEdge,orig_xBinLowEdge, orig_xBinUpEdge))

        #consistency check whether the bin-edges are the same between old and rebinned histo
        if (orig_xBinLowEdge != xBinLowEdge) or (orig_xBinUpEdge !=  xBinUpperEdge) :
            print "ERROR : The X bin edges of original and rebinned 2D histo are not the same!"
            print "Bin edges are original, rebinned = ", orig_xBinLowEdge, orig_xBinUpEdge, " : ", xBinLowEdge, xBinUpperEdge 
            sys.exit()

        for yBin in range(1,rebinnedHisto.GetNbinsY()+2):
            #Get bin-edges Y-axis
            yBinLowEdge=rebinnedHisto.GetYaxis().GetBinLowEdge(yBin)
            yBinUpperEdge=rebinnedHisto.GetYaxis().GetBinLowEdge(yBin)+rebinnedHisto.GetYaxis().GetBinWidth(yBin)

            #Find corresponding bins in original histogram
            orig_yBinLow = originalHisto.GetYaxis().FindBin(yBinLowEdge+0.001)
            orig_yBinUp = originalHisto.GetYaxis().FindBin(yBinUpperEdge-0.001)

            orig_yBinLowEdge = originalHisto.GetYaxis().GetBinLowEdge(orig_yBinLow)
            orig_yBinUpEdge = originalHisto.GetYaxis().GetBinLowEdge(orig_yBinUp)+originalHisto.GetYaxis().GetBinWidth(orig_yBinUp)
            # is correct
            if debug:
                print('New Bins Y [{:f}-{:f}] \t Original Bins  Y [{:f}-{:f}]'.format(yBinLowEdge, yBinUpperEdge, orig_yBinLowEdge, orig_yBinUpEdge))

            #consistency check whether the bin-edges are the same between old and rebinned histo
            if (isclose(orig_yBinLowEdge, yBinLowEdge, rel_tol=1e-3) != True  or (isclose(orig_yBinUpEdge, yBinUpperEdge, rel_tol=1e-3)!= True)  ) :
                print "ERROR : The Y bin edges of original and rebinned 2D histo are not the same!"
                print "Bin edges are original, rebinned = ", orig_yBinLowEdge,orig_yBinUpEdge, " : " , yBinLowEdge,yBinUpperEdge
                sys.exit()
            orig_error = ROOT.double(0)
            orig_integral = originalHisto.IntegralAndError(orig_xBinLow, orig_xBinUp, orig_yBinLow, orig_yBinUp, orig_error) 
            rebinnedHisto.SetBinContent(xBin,yBin,orig_integral/(abs(yBinUpperEdge-yBinLowEdge)*abs(xBinUpperEdge-xBinLowEdge)))
            rebinnedHisto.SetBinError(xBin,yBin,orig_error/(abs(yBinUpperEdge-yBinLowEdge)*abs(xBinUpperEdge-xBinLowEdge)))
            #rebinnedHisto.SetBinContent(xBin,yBin,orig_integral/(yBinUpperEdge-yBinLowEdge))
            #rebinnedHisto.SetBinError(xBin,yBin,orig_error/(yBinUpperEdge-yBinLowEdge))
            #rebinnedHisto.SetBinContent(xBin,yBin,orig_integral)
            #rebinnedHisto.SetBinError(xBin,yBin,orig_error)
    #rebinnedHisto.Draw("surf")
    #input()
    return rebinnedHisto

def DrawErrorBand(graph):
    isErrorBand = graph.GetErrorYhigh(0) != -1 and graph.GetErrorYlow(0) != -1
    npoints     = graph.GetN()
 
    if not isErrorBand:
        graph.Draw("l same")
        return
 
    # Declare individual TGraph objects used in drawing error band
    central, min, max = ROOT.TGraph(), ROOT.TGraph(), ROOT.TGraph()
    shapes = []
    for i in range((npoints-1)*4):
        shapes.append(ROOT.TGraph())
 
    # Set ownership of TGraph objects
    ROOT.SetOwnership(central, False)
    ROOT.SetOwnership(    min, False)
    ROOT.SetOwnership(    max, False)
    for shape in shapes:
        ROOT.SetOwnership(shape, False)
 
    # Get data points from TGraphAsymmErrors
    x, y, ymin, ymax = [], [], [], []
    for i in range(npoints):
        tmpX, tmpY = ROOT.Double(0), ROOT.Double(0)
        graph.GetPoint(i, tmpX, tmpY)
        x.append(tmpX)
        y.append(tmpY)
        ymin.append(tmpY - graph.GetErrorYlow(i))
        ymax.append(tmpY + graph.GetErrorYhigh(i))
 
    # Fill central, min and max graphs
    for i in range(npoints):
        central.SetPoint(i, x[i], y[i])
    min.SetPoint(i, x[i], ymin[i])
    max.SetPoint(i, x[i], ymax[i])
 
    # Fill shapes which will be shaded to create the error band
    for i in range(npoints-1):
        for version in range(4):
            shapes[i+(npoints-1)*version].SetPoint((version+0)%4, x[i],   ymax[i])
            shapes[i+(npoints-1)*version].SetPoint((version+1)%4, x[i+1], ymax[i+1])
            shapes[i+(npoints-1)*version].SetPoint((version+2)%4, x[i+1], ymin[i+1])
            shapes[i+(npoints-1)*version].SetPoint((version+3)%4, x[i],   ymin[i])
 
    # Set attributes to those of input graph
    central.SetLineColor(graph.GetLineColor())
    central.SetLineStyle(graph.GetLineStyle())
    central.SetLineWidth(graph.GetLineWidth())
    min.SetLineColor(graph.GetLineColor())
    min.SetLineStyle(graph.GetLineStyle())
    max.SetLineColor(graph.GetLineColor())
    max.SetLineStyle(graph.GetLineStyle())
    for shape in shapes:
        shape.SetFillColor(graph.GetFillColor())
        shape.SetFillStyle(graph.GetFillStyle())
 
    # Draw
    for shape in shapes:
        shape.Draw("f same")
    min.Draw("p same")
    max.Draw("p same")
    central.Draw("p same")
    ROOT.gPad.RedrawAxis()

def GetVariationEntryInFile(key):
    #Dictionary for variations in file. Number indicates entry in the Top Drawer file
    variationDict = {
    "1.000" : 1,
    "0.167" : 2,
    "0.250"  : 3,
    "0.500" : 4,
    "0.125" : 5,
    "2.000" : 6,
    "4.000" : 7,
    "6.000" : 8,
    "8.000" : 9
    }

    return variationDict.get(key)

def titleAndLabels(variable):
    
    titleDict = {
        "pt top" : "Transverse Momentum of Top Quark",
        "pt H" : "Transverse Momentum of Higgs Boson",
        "y top" : "Rapidity of Top Quark",
        "y H" : "Rapidity of Higgs Boson",
        "tt inv" : "Invariant Mass of t#bar{t} System",
        "inv mass" : "Invariant Mass of t#bar{t}H System",
        "total rate" : "Total t#bar{t}H Cross-Section",
    }

    xLabelDict = {
        "pt top" : "p^{Top}_{T} [GeV]",
        "pt H" : "p^{H}_{T} [GeV]",
        "y top" : "y^{Top}",
        "y H" : "y^{H}",
        "tt inv" : "M_{t#bar{t}} [GeV]",
        "inv mass" : "M_{t#bar{t}H} [GeV]",
        "total rate" : "#mu",
    }
    yLabelDict = {
        "pt top"    : "#frac{d#sigma}{dp^{Top}_{T}} [fb]",
        "pt H"      : "#frac{d#sigma}{dp^{H}_{T}} [fb]",
        "y top"     : "#frac{d#sigma}{dy^{Top}} [fb]",
        "y H"       : "#frac{d#sigma}{dy^{H}} [fb]",
        "tt inv"    : "#frac{d#sigma}{dM_{t#bar{t}}} [fb]",
        "inv mass"  : "#frac{d#sigma}{dM_{t#bar{t}}} [fb]",
        "total rate" : "#sigma_{t#bar{t}H} [fb]",
    }


    title = titleDict.get(variable) 
    xLabel = xLabelDict.get(variable) 
    yLabel = yLabelDict.get(variable) 
    return title, xLabel, yLabel
    
def getYBins(variable):
    yBinsTotalRate = []
    
    yBinsPT = array.array('d',np.arange(0., 300., 30.))
    yBinsPT = np.append(yBinsPT,array.array('d',[360.,420.,480.,540.,640.]))

    #yBinsPT= []

    yBinsPTHiggs = array.array('d',np.arange(0, 300., 30.))
    yBinsPTHiggs = np.append(yBinsPTHiggs,array.array('d',[360.,420.,480.,540.,600.,700.]))

    #yBinsPTHiggs = []
    yBinsRapidity = array.array('d',np.arange(-3.0, 3.4, 0.4))
    #yBinsRapidity = []

    yBinsTTinv = array.array('d',np.arange(340., 800., 20.))
    yBinsTTinv = np.append(yBinsTTinv,array.array('d',[800.,900.,1000.,1200.]))
    #yBinsTTinv = []

    yBinsTTHinv = array.array('d',np.arange(460., 880., 40.))
    yBinsTTHinv = np.append(yBinsTTHinv,array.array('d',[880.,980.,1100.,1300.,1600.,1900.]))
    #yBinsTTHinv = []
    yBinsDict = {
        "pt top" : yBinsPT,
        "pt H"   : yBinsPTHiggs,
        "y top"  : yBinsRapidity,
        "y H"    : yBinsRapidity,
        "tt inv" : yBinsTTinv,
        "inv mass" : yBinsTTHinv,
        "total rate" : yBinsTotalRate
    }

    return yBinsDict.get(variable) 

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def printCrossSections(inputPoleMass, inputRunningMass):

    #errors contains the relvative ucnertainties up/down
    errors, errorsUp, errorsDown = extractUncertaintiesFromHistos(inputPoleMass)
    errorsRunning, errorsUpRunning, errorsDownRunning = extractUncertaintiesFromHistos(inputRunningMass)

    print("Bin Edges \t Pole Mass [fb/bin] \t Running Mass [fb/bin]")
    nBins = len(errors)

    for xBin in range(nBins):
        lowBinEdge  = inputPoleMass[0].GetBinLowEdge(xBin+1)
        highBinEdge = inputPoleMass[0].GetBinLowEdge(xBin+1)+inputPoleMass[0].GetBinWidth(xBin+1)
        binValue    = inputPoleMass[0].GetBinContent(xBin+1)
        binWidth    = inputPoleMass[0].GetBinWidth(xBin+1)

        lowBinEdgeRunning  = inputRunningMass[0].GetBinLowEdge(xBin+1)
        highBinEdgeRunning = inputRunningMass[0].GetBinLowEdge(xBin+1)+inputRunningMass[0].GetBinWidth(xBin+1)
        binValueRunning    = inputRunningMass[0].GetBinContent(xBin+1)
        binWidthRunning    = inputRunningMass[0].GetBinWidth(xBin+1)

        print(" $[{:.2f},{:.2f}]$ \t & \t  ${{{:.3f}}}^{{+{:.3f}(+{:.2f}\\%)}}_{{{:.3f}({:.2f}\\%)}}$ & \t ${{{:.3f}}}^{{+{:.3f}(+{:.2f}\\%)}}_{{{:.3f}({:.2f}\\%)}}$  \\\\ ".format(lowBinEdge,\
         highBinEdge, binValue, binValue*errors[xBin][0], 100.*errors[xBin][0], binValue*errors[xBin][1], 100.*errors[xBin][1],\
         binValueRunning, binValueRunning*errorsRunning[xBin][0], 100.*errorsRunning[xBin][0], binValueRunning*errorsRunning[xBin][1], 100.*errorsRunning[xBin][1] ))

    print("total Cross Section Pole Mass: {:f} fb".format(inputPoleMass[0].Integral("width")))
    print("total Cross Section Running Mass: {:f} fb".format(inputRunningMass[0].Integral("width")))
    return