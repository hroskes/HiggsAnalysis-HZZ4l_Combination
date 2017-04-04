#! /usr/bin/env python
import sys
import os
import re
import math
import subprocess
from scipy.special import erf
import ROOT
from array import array
from systematicsClass import systematicsClass
from inputReader import inputReader
from helperstuff.combinehelpers import getdatatree, gettemplate, discriminantnames
from helperstuff import constants
from helperstuff.enums import Analysis, Production

## ------------------------------------
##  card and workspace class
## ------------------------------------

class properties_datacardClass:

    def __init__(self, analysis, production):
        self.ID_4mu = 1
        self.ID_4e  = 2
        self.ID_2e2mu = 3
        self.isFSR = True
        self.analysis = Analysis(analysis)
        self.production = Production(production)
        self.sqrts = 13

    def loadIncludes(self):
        import include

    def getHistFunc(self,newname,hist,Dx,Dy,Dz):
        fh = ROOT.FastHisto3D_d(hist)
        obslist = ROOT.RooArgList(Dx,Dy,Dz)
        fhfcn = ROOT.FastHisto3DFunc_d(newname,"",obslist,fh)
        return fhfcn

    def getChannelName(self):
        channelName = ""
        if (self.channel == self.ID_4mu): channelName = "4mu"
        elif (self.channel == self.ID_4e): channelName = "4e"
        elif (self.channel == self.ID_2e2mu): channelName = "2e2mu"
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"
        return channelName

    # cs x br function
    def makeXsBrFunction(self,signalProc,rrvMH):

        procName = "ggH"
        if(signalProc == 0): procName = "ggH" #dummy, when you sum up all the 5 chans
        if(signalProc == 1): procName = "ggH"
        if(signalProc == 2): procName = "qqH"
        if(signalProc == 3): procName = "WH"
        if(signalProc == 4): procName = "ZH"
        if(signalProc == 5): procName = "ttH"

        channelName = self.getChannelName()

        myCSWrhf = ROOT.HiggsCSandWidth()

        histXsBr = ROOT.TH1F("hsmxsbr_{0}_{1}".format(procName,channelName),"", 8905, 109.55, 1000.05)

        for i in range(1,8906):

            mHVal = histXsBr.GetBinCenter(i)
            BR = 0.0
            if (self.channel == self.ID_2e2mu):
                BR = myCSWrhf.HiggsBR(13,mHVal)
            else:
                BR = myCSWrhf.HiggsBR(12,mHVal)

            if (signalProc == 3 or signalProc == 4 or signalProc == 5):
                #overwrite BR if VH,ttH sample
                #these samples have inclusive Z decay
                BR = myCSWrhf.HiggsBR(11,mHVal)

            if (signalProc==0):
                totXs=0
                for ch in range(1,6):
                    totXs+=myCSWrhf.HiggsCS(ch, mHVal, self.sqrts)
                histXsBr.SetBinContent(i, totXs * BR)
            else:
                histXsBr.SetBinContent(i, myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts) * BR)

            #print '\nmakeXsBrFunction : procName=',procName,'   signalProc=',signalProc,'  mH (input)=',rrvMH.getVal(),
            #print '   CS=',myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts),'   BR=',BR

        rdhname = "rdhXsBr_{0}_{1}_{2}".format(procName,self.channel,self.sqrts)
        rdhXsBr = ROOT.RooDataHist(rdhname,rdhname, ROOT.RooArgList(rrvMH), histXsBr)

        return rdhXsBr

    # return trueVar if testStatement else return falseVar
    def getVariable(self,trueVar,falseVar,testStatement):

        if (testStatement):
            return trueVar
        else:
            return falseVar

    # main datacard and workspace function
    def makeCardsWorkspaces(self, theMH, theOutputDir, theInputs, theOptions):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        DEBUG = False
        self.mH = 125
        self.SMDsigCut = 1.
        self.SMDbkgCut = 1.
        self.lumi = self.production.dataluminosity
        self.sqrts = 13
        self.channel = theInputs['decayChannel']
        self.bkgMorph = theInputs['useCMS_zz4l_zjet']
        self.outputDir = theOutputDir
        self.templateDir = theOptions.templateDir
        self.dataAppendDir = theOptions.dataDirAppend
        self.sigmaVVai = theOptions.sigmaVVaiVal

        self.ggH_chan = theInputs['ggH']
        self.qqH_chan = theInputs['qqH']
        self.WH_chan = theInputs['WH']
        self.ZH_chan = theInputs['ZH']
        self.ttH_chan = theInputs['ttH']
        self.qqZZ_chan = theInputs['qqZZ']
        self.ggZZ_chan = theInputs['ggZZ']
        self.zjets_chan = theInputs['zjets']

        ## ---------------- SET PLOTTING STYLE ---------------- ##
        ROOT.setTDRStyle(True)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPadLeftMargin(0.16)

        ## ---------------- VARIABLES FOR LATER --------------- ##
        self.bUseCBnoConvolution = False
        ForXSxBR = False

        myCSW = ROOT.HiggsCSandWidth()

        ## ----------------- WIDTH AND RANGES ----------------- ##
        self.widthHVal =  myCSW.HiggsWidth(0,self.mH)
        if(self.widthHVal < 0.12):
            self.bUseCBnoConvolution = True
        self.isHighMass = False
        if self.mH >= 390:
            if theInputs['useHighMassReweightedShapes']:
                self.isHighMass = True
            else: print "useHighMassReweightedShapes set to FALSE, using non-reweighted shapes!"


        print "width: ",self.widthHVal

        self.windowVal = max( self.widthHVal, 1.0)
        lowside = 100.0
        highside = 1000.0

        if (self.mH >= 275):
            lowside = 180.0
            highside = 650.0
        if (self.mH >= 350):
            lowside = 200.0
            highside = 900.0
        if (self.mH >= 500):
            lowside = 250.0
            highside = 1000.0
        if (self.mH >= 700):
            lowside = 350.0
            highside = 1400.0

        self.low_M = max( (self.mH - 20.*self.windowVal), lowside)
        self.high_M = min( (self.mH + 15.*self.windowVal), highside)

        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"



        ## ------------------------- SYSTEMATICS CLASSES ----------------------------- ##

        systematics = systematicsClass( self.mH, False, self.isFSR, theInputs)
        systematics_forXSxBR = systematicsClass( self.mH, True, self.isFSR,theInputs)

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##

        bins = 1000
        if(self.bUseCBnoConvolution): bins = 200

        CMS_zz4l_mass_name = "ZZMass"

        CMS_zz4l_mass = ROOT.RooRealVar(CMS_zz4l_mass_name,CMS_zz4l_mass_name,self.low_M,self.high_M)
        CMS_zz4l_mass.setBins(bins)

        x_name = "CMS_zz4l_fai1"

        x = ROOT.RooRealVar(x_name,x_name,-1.,1.)
        x.setBins(bins)

        alpha_name = "CMS_zz4l_alpha"
        alpha_zz4l = ROOT.RooRealVar(alpha_name,alpha_name,0.,-1.,1.)
        alpha_zz4l.setBins(bins)


        D1Name, D2Name, D3Name = discriminantnames(self.analysis)

        self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.production.year),"LUMI_{0:.0f}".format(self.production.year),self.lumi)
        self.LUMI.setConstant(True)

        self.MH = ROOT.RooRealVar("MH","MH",self.mH)
        self.MH.setConstant(True)

        self.R = ROOT.RooRealVar("R","R",1.,0.,400.)
        self.RF = ROOT.RooRealVar("RF","RF",1.,0.,400.)
        self.RV = ROOT.RooRealVar("RV","RV",1.,0.,400.)
        Rsqrts_name = "R_{0:.0f}TeV".format(self.sqrts)
        RFsqrts_name = "RF_{0:.0f}TeV".format(self.sqrts)
        RVsqrts_name = "RV_{0:.0f}TeV".format(self.sqrts)
        self.Rsqrts = ROOT.RooRealVar(Rsqrts_name,Rsqrts_name,1.,0.,400.)
        self.RFsqrts = ROOT.RooRealVar(RFsqrts_name,RFsqrts_name,1.,0.,400.)
        self.RVsqrts = ROOT.RooRealVar(RVsqrts_name,RVsqrts_name,1.,0.,400.)
        self.muF = ROOT.RooFormulaVar("muF_{0:.0f}TeV".format(self.sqrts),"@0*@1*@2*@3",ROOT.RooArgList(self.R,self.Rsqrts,self.RF,self.RFsqrts))
        self.muV = ROOT.RooFormulaVar("muV_{0:.0f}TeV".format(self.sqrts),"@0*@1*@2*@3",ROOT.RooArgList(self.R,self.Rsqrts,self.RV,self.RVsqrts))

        self.sigmaVVaiVal = dict()
        for key, value in self.sigmaVVai.iteritems():
           rrvname = "{0}_{1:.0f}TeV".format(key,self.sqrts)
           rrv = ROOT.RooConstVar(rrvname,rrvname,value)
           self.sigmaVVaiVal[key] = rrv
        for key, value in self.sigmaVVaiVal.iteritems():
           print "{} = {}".format(key,value.getVal())

	# n2, alpha2 are right side parameters of DoubleCB
	# n, alpha are left side parameters of DoubleCB

        n_CB_d = 0.0
        alpha_CB_d = 0.0
        n2_CB_d = 0.0
        alpha2_CB_d = 0.0
        mean_CB_d = 0.0
        sigma_CB_d = 0.0
        mean_BW_d = self.mH
        gamma_BW_d = 0.0

        rdhXsBrFuncV_1 = self.makeXsBrFunction(1,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ggH",self.channel,self.production.year)
        rhfXsBrFuncV_1 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_1, 1)

        rdhXsBrFuncV_2 = self.makeXsBrFunction(2,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("VBF",self.channel,self.production.year)
        rhfXsBrFuncV_2 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_2, 1)

        rdhXsBrFuncV_3 = self.makeXsBrFunction(3,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("WH",self.channel,self.production.year)
        rhfXsBrFuncV_3 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_3, 1)

        rdhXsBrFuncV_4 = self.makeXsBrFunction(4,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ZH",self.channel,self.production.year)
        rhfXsBrFuncV_4 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_4, 1)

        rdhXsBrFuncV_5 = self.makeXsBrFunction(5,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ttH",self.channel,self.production.year)
        rhfXsBrFuncV_5 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_5, 1)


        ## -------- Variable Definitions -------- ##
        name = "CMS_zz4l_mean_e_sig"
        CMS_zz4l_mean_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_mean_e_err_{0}_{1:.0f}".format(self.channel,self.production.year)
        CMS_zz4l_mean_e_err = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_err",float(theInputs['CMS_zz4l_mean_e_sig']),-0.99,0.99)
        name = "CMS_zz4l_sigma_e_sig"
        CMS_zz4l_sigma_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_e_sig",3.0,0.0,30.0)
        name = "CMS_zz4l_mean_m_sig"
        CMS_zz4l_mean_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_mean_m_err_{0}_{1:.0f}".format(self.channel,self.production.year)
        CMS_zz4l_mean_m_err = ROOT.RooRealVar(name,"CMS_zz4l_mean_m_err",float(theInputs['CMS_zz4l_mean_m_sig']),-0.99,0.99)
        name = "CMS_zz4l_sigma_m_sig"
        CMS_zz4l_sigma_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)


        name = "CMS_zz4l_alpha2_{0}_{1:.0f}".format(self.channel,self.production.year)
        CMS_zz4l_alpha2 = ROOT.RooRealVar(name,"CMS_zz4l_alpha2",1.,-10.,10.)
        name = "CMS_zz4l_n2_sig_{0}_{1:.0f}".format(self.channel,self.production.year)
        CMS_zz4l_n2 = ROOT.RooRealVar(name,"CMS_zz4l_n2",2.,-10.,10.)
        name = "CMS_zz4l_alpha_{0}_{1:.0f}".format(self.channel,self.production.year)
        CMS_zz4l_alpha = ROOT.RooRealVar(name,"CMS_zz4l_alpha",1.,-10.,10.)
        name = "CMS_zz4l_n_sig_{0}_{1:.0f}".format(self.channel,self.production.year)
        CMS_zz4l_n = ROOT.RooRealVar(name,"CMS_zz4l_n",2.,-10.,10.)
        name = "CMS_zz4l_mean_BW_{0}_{1:.0f}".format(self.channel,self.production.year)
        CMS_zz4l_mean_BW = ROOT.RooRealVar(name,"CMS_zz4l_mean_BW",self.mH,self.low_M,self.high_M)
        name = "interf_ggH"
        #name = "CMS_zz4l_gamma_sig_{0}_{1:.0f}".format(self.channel,self.production.year)
        CMS_zz4l_gamma = ROOT.RooRealVar(name,"CMS_zz4l_gamma",10.,0.001,1000.)
        name = "CMS_zz4l_widthScale_{0}_{1:.0f}".format(self.channel,self.production.year)
        CMS_zz4l_widthScale = ROOT.RooRealVar(name,"CMS_zz4l_widthScale",1.0)

        one = ROOT.RooRealVar("one","one",1.0)
        one.setConstant(True)

        CMS_zz4l_mean_BW.setVal( mean_BW_d )
        CMS_zz4l_gamma.setVal(0)
        CMS_zz4l_mean_e_sig.setVal(0)
        CMS_zz4l_mean_e_err.setConstant(True)
        CMS_zz4l_sigma_e_sig.setVal(0)
        CMS_zz4l_mean_m_sig.setVal(0)
        CMS_zz4l_mean_m_err.setConstant(True)
        CMS_zz4l_sigma_m_sig.setVal(0)
        CMS_zz4l_alpha.setVal(0)
        CMS_zz4l_n.setVal(0)
        CMS_zz4l_alpha2.setVal(0)
        CMS_zz4l_n2.setVal(0)

        CMS_zz4l_widthScale.setConstant(True)
        CMS_zz4l_mean_BW.setConstant(True)

        print "mean_BW ", CMS_zz4l_mean_BW.getVal()
        print "gamma_BW ", CMS_zz4l_gamma.getVal()
        print "mean_e_sig ", CMS_zz4l_mean_e_sig.getVal()
        print "mean_e_err ", CMS_zz4l_mean_e_err.getVal()
        print "sigma_e ", CMS_zz4l_sigma_e_sig.getVal()
        print "mean_m_sig ",CMS_zz4l_mean_m_sig.getVal()
        print "mean_m_err ", CMS_zz4l_mean_m_err.getVal()
        print "sigma_m ", CMS_zz4l_sigma_m_sig.getVal()
        print "alpha ", CMS_zz4l_alpha.getVal()
        print "n ", CMS_zz4l_n.getVal()
        print "alpha2 ", CMS_zz4l_alpha2.getVal()
        print "n2 ", CMS_zz4l_n2.getVal()




        ## -------------------- RooFormulaVar's -------------------- ##
        rfv_n_CB = ROOT.RooFormulaVar()
        rfv_alpha_CB = ROOT.RooFormulaVar()
        rfv_n2_CB = ROOT.RooFormulaVar()
        rfv_alpha2_CB = ROOT.RooFormulaVar()
        rfv_mean_CB = ROOT.RooFormulaVar()
        rfv_sigma_CB = ROOT.RooFormulaVar()

        name = "CMS_zz4l_n_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.production.year)
        if self.isHighMass : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
        else : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))

        name = "CMS_zz4l_alpha_{0:.0f}_centralValue".format(self.channel)
        if self.isHighMass : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape_HM'], ROOT.RooArgList(self.MH))
        else : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape'], ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_n2_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.production.year)
        #if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
        #else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
        if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")",ROOT.RooArgList(self.MH))
        else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")",ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_alpha2_{0:.0f}_centralValue".format(self.channel)
        if self.isHighMass : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape_HM'], ROOT.RooArgList(self.MH))
        else : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape'], ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_mean_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.production.year)


        if (self.channel == self.ID_4mu) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_m_err))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_m_err))
        elif (self.channel == self.ID_4e) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig,CMS_zz4l_mean_e_err))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig,CMS_zz4l_mean_e_err))
        elif (self.channel == self.ID_2e2mu) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+ (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig,CMS_zz4l_mean_m_err,CMS_zz4l_mean_e_err))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+ (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig,CMS_zz4l_mean_m_err,CMS_zz4l_mean_e_err))


        name = "CMS_zz4l_sigma_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.production.year)

        if (self.channel == self.ID_4mu) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*TMath::Sqrt((1+@1)*(1+@2))", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*TMath::Sqrt((1+@1)*(1+@2))", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))

        name = "CMS_zz4l_gamma_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.production.year)
        rfv_gamma_BW = ROOT.RooFormulaVar(name,"("+theInputs['gamma_BW_shape_HM']+")"+"*(1+@1*0.05)",ROOT.RooArgList(self.MH,CMS_zz4l_gamma))

        if (DEBUG): print " DEBUG *********  ", theInputs['sigma_CB_shape']

        print "n_CB ", rfv_n_CB.getVal()
        print "alpha_CB ", rfv_alpha_CB.getVal()
        print "n2_CB ", rfv_n2_CB.getVal()
        print "alpha2_CB ", rfv_alpha2_CB.getVal()
        print "mean_CB ", rfv_mean_CB.getVal()
        print "sigma_CB ", rfv_sigma_CB.getVal()
        print "gamma_BW ", rfv_gamma_BW.getVal()

        CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_{0:.0f}_{1:.0f}".format(self.channel,self.production.year),"@0+@1", ROOT.RooArgList(rfv_mean_CB, self.MH))

        print "mean_sig_NoConv ", CMS_zz4l_mean_sig_NoConv.getVal()



        ## --------------------- SHAPE FUNCTIONS ---------------------- ##

        signalCB_ggH = ROOT.RooDoubleCB("signalCB_ggH","signalCB_ggH",CMS_zz4l_mass, self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB, self.bUseCBnoConvolution) , rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ggH = ROOT.RooRelBWUFParam("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ggH =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH,signalCB_ggH, 2)
        #High mass pdf
        signalBW_ggH_HM = ROOT.RooRelBWHighMass("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ggH_HM =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH_HM,signalCB_ggH, 2)


        signalCB_VBF = ROOT.RooDoubleCB("signalCB_VBF","signalCB_VBF",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_VBF = ROOT.RooRelBWUFParam("signalBW_VBF", "signalBW_VBF",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_VBF = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",CMS_zz4l_mass,signalBW_VBF,signalCB_VBF, 2)
        #High mass pdf
        signalBW_VBF_HM = ROOT.RooRelBWHighMass("signalBW_VBF", "signalBW_VBF",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_VBF_HM = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",CMS_zz4l_mass,signalBW_VBF_HM,signalCB_VBF, 2)


        signalCB_WH = ROOT.RooDoubleCB("signalCB_WH","signalCB_WH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_WH = ROOT.RooRelBWUFParam("signalBW_WH", "signalBW_WH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_WH = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",CMS_zz4l_mass,signalBW_WH,signalCB_WH, 2)
        #High mass pdf
        signalBW_WH_HM = ROOT.RooRelBWHighMass("signalBW_WH", "signalBW_WH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_WH_HM = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",CMS_zz4l_mass,signalBW_WH_HM,signalCB_WH, 2)


        signalCB_ZH = ROOT.RooDoubleCB("signalCB_ZH","signalCB_ZH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ZH = ROOT.RooRelBWUFParam("signalBW_ZH", "signalBW_ZH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ZH = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",CMS_zz4l_mass,signalBW_ZH,signalCB_ZH, 2)
        #High mass pdf
        signalBW_ZH_HM = ROOT.RooRelBWHighMass("signalBW_ZH", "signalBW_ZH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ZH_HM = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",CMS_zz4l_mass,signalBW_ZH_HM,signalCB_ZH, 2)


        signalCB_ttH = ROOT.RooDoubleCB("signalCB_ttH","signalCB_ttH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ttH = ROOT.RooRelBWUFParam("signalBW_ttH", "signalBW_ttH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ttH = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",CMS_zz4l_mass,signalBW_ttH,signalCB_ttH, 2)
        #High mass pdf
        signalBW_ttH_HM = ROOT.RooRelBWHighMass("signalBW_ttH", "signalBW_ttH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ttH_HM = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",CMS_zz4l_mass,signalBW_ttH_HM,signalCB_ttH, 2)


        ## Buffer fraction for cyclical behavior
        sig_ggH.setBufferFraction(0.2)
        sig_VBF.setBufferFraction(0.2)
        sig_WH.setBufferFraction(0.2)
        sig_ZH.setBufferFraction(0.2)
        sig_ttH.setBufferFraction(0.2)

        sig_ggH_HM.setBufferFraction(0.2)
        sig_VBF_HM.setBufferFraction(0.2)
        sig_WH_HM.setBufferFraction(0.2)
        sig_ZH_HM.setBufferFraction(0.2)
        sig_ttH_HM.setBufferFraction(0.2)

        ## -------------------- 1D SIGNAL SHAPES FOR PROPERTIES ------------------------- ##

        print '1D signal shapes for Properties'

        channelName = ""
        if (self.channel == self.ID_4mu): channelName = "4mu"
        elif (self.channel == self.ID_4e): channelName = "4e"
        elif (self.channel == self.ID_2e2mu): channelName = "2e2mu"
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"

        Sig_T = [gettemplate(self.analysis, self.production, s, channelName) for s in self.analysis.signalsamples()]
        for i, t in enumerate(Sig_T):
            t.SetName("T_ZZ_{:.0f}_{}_3D_{}".format(self.production.year,self.appendName,i))

        Sig_T_ScaleResUp = [gettemplate(self.analysis, self.production, s, channelName, "ScaleResUp") for s in self.analysis.signalsamples()]
        for i, t in enumerate(Sig_T_ScaleResUp):
            t.SetName("T_ZZ_{:.0f}_{}_3D_{}_ScaleResUp".format(self.production.year,self.appendName,i))

        Sig_T_ScaleResDown = [gettemplate(self.analysis, self.production, s, channelName, "ScaleResDown") for s in self.analysis.signalsamples()]
        for i, t in enumerate(Sig_T_ScaleResDown):
            t.SetName("T_ZZ_{:.0f}_{}_3D_{}_ScaleResDown".format(self.production.year,self.appendName,i))

        dBinsX = Sig_T[0].GetXaxis().GetNbins()
        print "X bins: ",dBinsX
        dLowX = Sig_T[0].GetXaxis().GetXmin()
        dHighX = Sig_T[0].GetXaxis().GetXmax()

        dBinsY = Sig_T[0].GetYaxis().GetNbins()
        print "Y bins: ",dBinsY
        dLowY = Sig_T[0].GetYaxis().GetXmin()
        dHighY = Sig_T[0].GetYaxis().GetXmax()

        dBinsZ = Sig_T[0].GetZaxis().GetNbins()
        print "Z bins: ",dBinsZ
        dLowZ = Sig_T[0].GetZaxis().GetXmin()
        dHighZ = Sig_T[0].GetZaxis().GetXmax()

        D1 = ROOT.RooRealVar(D1Name,D1Name,dLowX,dHighX)
        D2 = ROOT.RooRealVar(D2Name,D2Name,dLowY,dHighY)
        D3 = ROOT.RooRealVar(D3Name,D3Name,dLowZ,dHighZ)
        D1.setBins(dBinsX)
        D2.setBins(dBinsY)
        D3.setBins(dBinsZ)

        SigHFcn = [self.getHistFunc("{}_histfunc".format(tpl.GetName()),tpl,D1,D2,D3) for tpl in Sig_T]
        SigHFcn_SRUp = [self.getHistFunc("{}_histfunc".format(tpl.GetName()),tpl,D1,D2,D3) for tpl in Sig_T_ScaleResUp]
        SigHFcn_SRDown = [self.getHistFunc("{}_histfunc".format(tpl.GetName()),tpl,D1,D2,D3) for tpl in Sig_T_ScaleResDown]

        ggHpdfName = "ggH_RooSpinZeroPdf_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ggHpdf = ROOT.HZZ4L_RooSpinZeroPdf_1D_fast(
           ggHpdfName,ggHpdfName,
           x,
           ROOT.RooArgList(D1,D2,D3),
           ROOT.RooArgList(
              SigHFcn[0],SigHFcn[1],SigHFcn[2]
           )
        )

        ggHpdfName_syst1Up = "ggH_RooSpinZeroPdf_ScaleResUp_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ggHpdf_syst1Up = ROOT.HZZ4L_RooSpinZeroPdf_1D_fast(
           ggHpdfName_syst1Up,ggHpdfName_syst1Up,
           x,
           ROOT.RooArgList(D1,D2,D3),
           ROOT.RooArgList(
              SigHFcn_SRUp[0],SigHFcn_SRUp[1],SigHFcn_SRUp[2]
           )
        )

        ggHpdfName_syst1Down = "ggH_RooSpinZeroPdf_ScaleResDown_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ggHpdf_syst1Down = ROOT.HZZ4L_RooSpinZeroPdf_1D_fast(
           ggHpdfName_syst1Down,ggHpdfName_syst1Down,
           x,
           ROOT.RooArgList(D1,D2,D3),
           ROOT.RooArgList(
              SigHFcn_SRDown[0],SigHFcn_SRDown[1],SigHFcn_SRDown[2]
           )
        )


        ## ------------------ END 2D SIGNAL SHAPES FOR PROPERTIES ------------------------ ##


        ## ------------------ 2D BACKGROUND SHAPES FOR PROPERTIES ------------------- ##

        qqZZTemplate = gettemplate(self.analysis, self.production, "qqZZ", channelName)

        TemplateName = "qqZZTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        qqZZTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),qqZZTemplate)
        PdfName = "qqZZ_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        qqZZTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),qqZZTempDataHist)

        ggZZTemplate = gettemplate(self.analysis, self.production, "ggZZ", channelName)

        TemplateName = "ggZZTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ggZZTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),ggZZTemplate)
        PdfName = "ggZZ_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ggZZTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),ggZZTempDataHist)


        ZjetsTemplate = gettemplate(self.analysis, self.production, "ZX", channelName)
        TemplateName = "ZjetsTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ZjetsTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),ZjetsTemplate)
        PdfName = "Zjets_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ZjetsTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),ZjetsTempDataHist)

        ZjetsTemplateDown = gettemplate(self.analysis, self.production, "ZX", channelName, "ZXDown")
        TemplateName = "ZjetsTempDownDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ZjetsTempDataHistDown = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),ZjetsTemplateDown)
        PdfName = "Zjets_TemplateDownPdf_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ZjetsTemplatePdfDown = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),ZjetsTempDataHistDown)

        ZjetsTemplateUp = gettemplate(self.analysis, self.production, "ZX", channelName, "ZXUp")
        TemplateName = "ZjetsTempUpDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ZjetsTempDataHistUp = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),ZjetsTemplateUp)
        PdfName = "Zjets_TemplateUpPdf_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ZjetsTemplatePdfUp = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),ZjetsTempDataHistUp)

        funcList_zjets = ROOT.RooArgList()
        morphBkgVarName =  "CMS_zz4l_smd_zjets_bkg_{}".format(self.getChannelName())
        alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        morphVarListBkg = ROOT.RooArgList()

        if(self.bkgMorph):
            funcList_zjets.add(ZjetsTemplatePdf)
            funcList_zjets.add(ZjetsTemplatePdfUp)
            funcList_zjets.add(ZjetsTemplatePdfDown)
            alphaMorphBkg.setConstant(False)
            morphVarListBkg.add(alphaMorphBkg)
        else:
            funcList_zjets.add(ZjetsTemplatePdf)
            alphaMorphBkg.setConstant(True)
        MorphName = "ZX_TemplateMorphPdf_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        ZjetsTemplateMorphPdf = ROOT.FastVerticalInterpHistPdf3D(MorphName,MorphName,D1,D2,D3,False,funcList_zjets,morphVarListBkg,1.0,1)


        ## ---------------- END 2D BACKGROUND SHAPES FOR PROPERTIES ----------------- ##

        ## ------------------- LUMI -------------------- ##

        rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)

        ## ----------------------- SIGNAL RATES ----------------------- ##

        CS_ggH = constants.SMXSggH #take from YR4
        CS_VBF = constants.SMXSVBF
        CS_WH = constants.SMXSWH
        CS_ZH = constants.SMXSZH
        CS_ttH = constants.SMXSttH

        CS_ggH_rrvname = "CSggHval_{0:.0f}".format(self.sqrts)
        CS_VBF_rrvname = "CSVBFval_{0:.0f}".format(self.sqrts)
        CS_WH_rrvname = "CSWHval_{0:.0f}".format(self.sqrts)
        CS_ZH_rrvname = "CSZHval_{0:.0f}".format(self.sqrts)
        CS_ttH_rrvname = "CSttHval_{0:.0f}".format(self.sqrts)
        CS_total_rrvname = "CStotalval_{0:.0f}".format(self.sqrts)
        CS_totalff_rrvname = "CStotalffval_{0:.0f}".format(self.sqrts)
        CS_totalVV_rrvname = "CStotalVVval_{0:.0f}".format(self.sqrts)
        CS_fracff_rrvname = "CSfracffval_{0:.0f}".format(self.sqrts)
        CS_fracVV_rrvname = "CSfracVVval_{0:.0f}".format(self.sqrts)
        self.CS_ggH_rrv = ROOT.RooConstVar(CS_ggH_rrvname,CS_ggH_rrvname,CS_ggH)
        self.CS_VBF_rrv = ROOT.RooConstVar(CS_VBF_rrvname,CS_VBF_rrvname,CS_VBF)
        self.CS_WH_rrv = ROOT.RooConstVar(CS_WH_rrvname,CS_WH_rrvname,CS_WH)
        self.CS_ZH_rrv = ROOT.RooConstVar(CS_ZH_rrvname,CS_ZH_rrvname,CS_ZH)
        self.CS_ttH_rrv = ROOT.RooConstVar(CS_ttH_rrvname,CS_ttH_rrvname,CS_ttH)
        self.CSfflist = ROOT.RooArgList(self.CS_ggH_rrv,self.CS_ttH_rrv)
        self.CSVVlist = ROOT.RooArgList(self.CS_VBF_rrv,self.CS_WH_rrv,self.CS_ZH_rrv)
        self.CStotalff = ROOT.RooFormulaVar(CS_totalff_rrvname,"@0+@1",self.CSfflist)
        self.CStotalVV = ROOT.RooFormulaVar(CS_totalVV_rrvname,"@0+@1+@2",self.CSVVlist)
        self.CStotal = ROOT.RooFormulaVar(CS_total_rrvname,"@0+@1",ROOT.RooArgList(self.CStotalff,self.CStotalVV))
        self.CSfracff = ROOT.RooFormulaVar(CS_fracff_rrvname,"@0/@1",ROOT.RooArgList(self.CStotalff,self.CStotal))
        self.CSfracVV = ROOT.RooFormulaVar(CS_fracVV_rrvname,"@0/@1",ROOT.RooArgList(self.CStotalVV,self.CStotal))

        T1_integralName = "normt1_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        T2_integralName = "normt2_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        T4_integralName = "normt4_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        T1_integral = ROOT.RooConstVar (T1_integralName,T1_integralName,Sig_T[0].Integral())
        T2_integral = ROOT.RooConstVar (T2_integralName,T2_integralName,Sig_T[1].Integral())
        T4_integral = ROOT.RooConstVar (T4_integralName,T4_integralName,Sig_T[2].Integral())
        print "T1 integral",T1_integral.getVal()
        print "T2 integral",T2_integral.getVal()
        print "T4 integral",T4_integral.getVal()
        r_fai_pures_norm_Name = "sig_PuresNorm_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        r_fai_realints_norm_Name = "sig_RealIntsNorm_{0:.0f}_{1:.0f}".format(self.channel,self.production.year)
        r_fai_pures_norm = ROOT.RooFormulaVar(r_fai_pures_norm_Name,"( (1-abs(@0))*@1+abs(@0)*@2 )/@1",RooArgList(x,T1_integral,T2_integral))
        r_fai_realints_norm = ROOT.RooFormulaVar(r_fai_realints_norm_Name,"( sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*@1 )/@2",RooArgList(x,T4_integral,T1_integral))

        self.r_fai_norm = None
        self.r_fai_norm_dec = None
        #rf_fai_norm_prod = None
        self.rv_fai_norm_prod = None
        self.rv_fai_pures_norm = None
        self.rv_fai_realints_norm = None
        if theOptions.newMu:
          self.r_fai_norm_dec_name = "sig_DecNormPar_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
          self.r_fai_norm_dec = ROOT.RooFormulaVar(self.r_fai_norm_dec_name,"TMath::Max((@0+@1)*(1-abs(@2)),0)",ROOT.RooArgList(r_fai_pures_norm,r_fai_realints_norm,alpha_zz4l))

          self.rv_fai_pures_norm_Name = "sig_VV_PuresNorm_{0:.0f}".format(self.sqrts)
          self.rv_fai_realints_norm_Name = "sig_VV_RealIntsNorm_{0:.0f}".format(self.sqrts)
          self.rv_fai_norm_prod_Name = "sig_VV_Total_{0:.0f}".format(self.sqrts)

          self.rv_fai_pures_norm = ROOT.RooFormulaVar(self.rv_fai_pures_norm_Name,"( (1-abs(@0))*@1+abs(@0)*@2 )/@1",ROOT.RooArgList(x,self.sigmaVVaiVal["T1"],self.sigmaVVaiVal["T2"]))
          self.rv_fai_realints_norm = ROOT.RooFormulaVar(self.rv_fai_realints_norm_Name,"( sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*@1 )/@2",ROOT.RooArgList(x,self.sigmaVVaiVal["T4"],self.sigmaVVaiVal["T1"]))
          self.rv_fai_norm_prod = ROOT.RooFormulaVar(self.rv_fai_norm_prod_Name,"TMath::Max((@0+@1),0)",ROOT.RooArgList(self.rv_fai_pures_norm,self.rv_fai_realints_norm))

          self.r_fai_norm = ROOT.RooFormulaVar("ggH_norm","@0*(@1*@2 + @3*@4*@5)",ROOT.RooArgList(self.r_fai_norm_dec, self.muF,self.CSfracff, self.muV,self.rv_fai_norm_prod,self.CSfracVV))
        else:
          self.r_fai_norm = ROOT.RooFormulaVar("ggH_norm","TMath::Max((@0+@1)*(1-abs(@2)),0)",ROOT.RooArgList(r_fai_pures_norm,r_fai_realints_norm,alpha_zz4l))


        ## --------------------------- DATASET --------------------------- ##

        data_obs_tree = getdatatree(channelName, self.production)
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_{0}".format(self.appendName)


        data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,D1,D2,D3))


        ## --------------------------- WORKSPACE -------------------------- ##

        endsInP5 = False
        tmpMH = self.mH
        if ( math.fabs(math.floor(tmpMH)-self.mH) > 0.001): endsInP5 = True
        if (DEBUG): print "ENDS IN P5  ",endsInP5
        name_Shape = ""
        name_ShapeWS = ""
        name_ShapeWSXSBR = ""

        name_Shape2 = "hzz4l_{0}S_{1:.0f}.txt".format(self.appendName,self.production.year)
        name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}.input.root".format(self.appendName,self.production.year)
        if (endsInP5):
           name_Shape       =       "{0}/HCG/{1:.1f}/{2:.0f}/{3}".format(self.outputDir,self.mH,self.production.year,name_Shape2)
           name_ShapeWS     =       "{0}/HCG/{1:.1f}/{2:.0f}/{3}".format(self.outputDir,self.mH,self.production.year,name_ShapeWS2)
           name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.1f}/{2:.0f}/{3}".format(self.outputDir,self.mH,self.production.year,name_ShapeWS2)
        else:
           name_Shape       =       "{0}/HCG/{1:.0f}/{2:.0f}/{3}".format(self.outputDir,self.mH,self.production.year,name_Shape2)
           name_ShapeWS     =       "{0}/HCG/{1:.0f}/{2:.0f}/{3}".format(self.outputDir,self.mH,self.production.year,name_ShapeWS2)
           name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/{2:.0f}/{3}".format(self.outputDir,self.mH,self.production.year,name_ShapeWS2)
        try:
            os.makedirs(os.path.dirname(name_Shape))
        except OSError:
            pass
        try:
            os.makedirs(os.path.dirname(name_ShapeWSXSBR))
        except OSError:
            pass
        if(DEBUG): print name_Shape,"  ",name_ShapeWS2

        w = ROOT.RooWorkspace("w","w")
        #w.importClassCode(ROOT.RooqqZZPdf_v2.Class(),True)
        #w.importClassCode(ROOT.RooggZZPdf_v2.Class(),True)
        #w.importClassCode(ROOT.FastHistoFunc_f.Class(),True)
        #w.importClassCode(ROOT.FastHisto2DFunc_f.Class(),True)
        w.importClassCode(ROOT.FastHisto3DFunc_f.Class(),True)
        w.importClassCode(ROOT.HZZ4L_RooSpinZeroPdf_1D_fast.Class(),True)
        w.importClassCode(ROOT.RooFormulaVar.Class(),True)
        if self.isHighMass :
            w.importClassCode(ROOT.RooRelBWHighMass.Class(),True)


        getattr(w,'import')(data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?

        if self.r_fai_norm is None:
           print "ERROR: self.r_fai_norm is None!"
           sys.exit()
        else:
           print "Importing {}".format(self.r_fai_norm.GetName())
           self.r_fai_norm.SetName("ggH_norm")
           getattr(w,'import')(self.r_fai_norm,ROOT.RooCmdArg()) ### Should this be renamed?
           #getattr(w,'import')(self.r_fai_norm, ROOT.RooFit.Rename("ggH_norm")) ### Should this be renamed?
           self.r_fai_norm.Print("v")

        if theOptions.newMu:
           print "Importing new mu parameterization variables"
           getattr(w,'import')(self.R, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.RF, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.RV, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.Rsqrts, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.RFsqrts, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.RVsqrts, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.muF, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.muV, ROOT.RooFit.RecycleConflictNodes())
           for key, value in self.sigmaVVaiVal.iteritems():
              print "\tImporting {}={}".format(key,value.GetName())
              getattr(w,'import')(value, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CS_ggH_rrv, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CS_VBF_rrv, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CS_WH_rrv, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CS_ZH_rrv, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CS_ttH_rrv, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CStotalff, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CStotalVV, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CStotal, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CSfracff, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.CSfracVV, ROOT.RooFit.RecycleConflictNodes())
           #getattr(w,'import')(self.r_fai_norm_dec,ROOT.RooCmdArg())
           getattr(w,'import')(self.rv_fai_pures_norm, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.rv_fai_realints_norm, ROOT.RooFit.RecycleConflictNodes())
           getattr(w,'import')(self.rv_fai_norm_prod, ROOT.RooFit.RecycleConflictNodes())

        ggHpdf.SetNameTitle("ggH","ggH")
        getattr(w,'import')(ggHpdf, ROOT.RooFit.RecycleConflictNodes())
        ggHpdf_syst1Up.SetNameTitle("ggH_Res{0}Up".format(self.appendName),"ggH_Res{0}Up".format(self.appendName))
        getattr(w,'import')(ggHpdf_syst1Up, ROOT.RooFit.RecycleConflictNodes())
        ggHpdf_syst1Down.SetNameTitle("ggH_Res{0}Down".format(self.appendName),"ggH_Res{0}Down".format(self.appendName))
        getattr(w,'import')(ggHpdf_syst1Down, ROOT.RooFit.RecycleConflictNodes())
        #ggHpdf_syst2Up.SetNameTitle("ggH_Scale{0}Up".format(self.appendName),"ggH_Scale{0}Up".format(self.appendName))
        #getattr(w,'import')(ggHpdf_syst2Up, ROOT.RooFit.RecycleConflictNodes())
        #ggHpdf_syst2Down.SetNameTitle("ggH_Scale{0}Down".format(self.appendName),"ggH_Scale{0}Down".format(self.appendName))
        #getattr(w,'import')(ggHpdf_syst2Down, ROOT.RooFit.RecycleConflictNodes())
        #getattr(w,'import')(Sig_T_1,"T_1")
        #getattr(w,'import')(Sig_T_2,"T_2")
        #getattr(w,'import')(Sig_T_4,"T_3")
        #getattr(w,'import')(Sig_T_1,"T_1_{0}_{1}".format(self.appendName,self.production.year))
        #getattr(w,'import')(Sig_T_2,"T_2_{0}_{1}".format(self.appendName,self.production.year))
        #getattr(w,'import')(Sig_T_4,"T_3_{0}_{1}".format(self.appendName,self.production.year))

        qqZZTemplatePdf.SetNameTitle("bkg_qqzz","bkg_qqzz")
        ggZZTemplatePdf.SetNameTitle("bkg_ggzz","bkg_ggzz")
        #ZjetsTemplatePdf.SetNameTitle("bkg_zjets","bkg_zjets")
        ZjetsTemplateMorphPdf.SetNameTitle("bkg_zjets","bkg_zjets")
        getattr(w,'import')(qqZZTemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(ggZZTemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(ZjetsTemplateMorphPdf, ROOT.RooFit.RecycleConflictNodes())
        #getattr(w,'import')(ZjetsTemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        #qqZZTemplateMorphPdf.SetNameTitle("bkg_qqzzMorph","bkg_qqzzMorph")
        #getattr(w,'import')(qqZZTemplateMorphPdf, ROOT.RooFit.RecycleConflictNodes())
        #testpdf.SetNameTitle("testpdf","testpdf")
        #getattr(w,'import')(testpdf, ROOT.RooFit.RecycleConflictNodes())


        w.writeToFile(name_ShapeWS)
        w.writeToFile(name_ShapeWSXSBR)

        ## --------------------------- DATACARDS -------------------------- ##

        ## If the channel is not declared in inputs, set rate = 0
        if not self.ggH_chan:  sigRate_ggH_Shape = 0
        if not self.qqH_chan:  sigRate_VBF_Shape = 0
        if not self.WH_chan:   sigRate_WH_Shape = 0
        if not self.ZH_chan:   sigRate_ZH_Shape = 0
        if not self.ttH_chan:  sigRate_ttH_Shape = 0

        if not self.qqZZ_chan:  bkgRate_qqzz_Shape = 0
        if not self.ggZZ_chan:  bkgRate_ggzz_Shape = 0
        if not self.zjets_chan: bkgRate_zjets_Shape = 0

        rates = {}
        rates['ggH'] = sigRate_ggH_Shape
        rates['qqH'] = sigRate_VBF_Shape
        rates['WH']  = sigRate_WH_Shape
        rates['ZH']  = sigRate_ZH_Shape
        rates['ttH'] = sigRate_ttH_Shape

        rates['qqZZ']  = 6
        rates['ggZZ']  = 6
        rates['zjets'] = 6
        rates['ttbar'] = 0
        rates['zbb']   = 0

        systematics.setSystematics(None, None, None)
        systematics_forXSxBR.setSystematics(None, None, None)

        ## Write Datacards
        fo = open( name_Shape, "wb")
        self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs.numEntries())

        systematics.WriteSystematics(fo, theInputs)
        systematics.WriteShapeSystematics(fo,theInputs)


        fo.close()



        ## forXSxBR

        fo = open( name_Shape, "wb" )

        self.WriteDatacard(fo, theInputs,name_ShapeWS2, rates, data_obs.numEntries())

        systematics_forXSxBR.WriteSystematics(fo, theInputs)
        systematics_forXSxBR.WriteShapeSystematics(fo,theInputs)
        fo.close()




    def WriteDatacard(self,file,theInputs,nameWS,theRates,obsEvents):

        numberSig = self.numberOfSigChan(theInputs)
        numberBg  = self.numberOfBgChan(theInputs)

        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")

        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS w:$PROCESS_$SYSTEMATIC\n".format(nameWS))
        file.write("------------\n")


        file.write("bin a{0} \n".format(self.channel))
        file.write("observation {0} \n".format(obsEvents))

        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.low_M,self.high_M))
        file.write("bin ")

        channelList=['ggH','qqH','WH','ZH','ttH','qqZZ','ggZZ','zjets','ttbar','zbb']

        channelName=['ggH','qqH','WH','ZH','ttH','bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']

        for chan in channelList:
            if theInputs[chan]:
                file.write("a{0} ".format(self.channel))
        file.write("\n")

        file.write("process ")

        i=0

        for chan in channelList:
            #print 'checking if ',chan,' is in the list of to-do'
            #print "{0} ".format(channelName[i])
            if theInputs[chan]:
                file.write("{0} ".format(channelName[i]))
                #print 'writing in card index=',i,'  chan=',chan
                #print "{0} ".format(channelName[i])
            i+=1


        file.write("\n")

        processLine = "process "

        for x in range(-numberSig+1,1):
            processLine += "{0} ".format(x)

        for y in range(1,numberBg+1):
            processLine += "{0} ".format(y)

        file.write(processLine)
        file.write("\n")

        file.write("rate ")
        for chan in channelList:
            if theInputs[chan]:
                file.write("{0:.4f} ".format(theRates[chan]))
        file.write("\n")
        file.write("------------\n")



    def numberOfSigChan(self,inputs):

        counter=0

        if inputs['ggH']: counter+=1
        if inputs['qqH']: counter+=1
        if inputs['WH']:  counter+=1
        if inputs['ZH']:  counter+=1
        if inputs['ttH']: counter+=1

        return counter

    def numberOfBgChan(self,inputs):

        counter=0

        if inputs['qqZZ']:  counter+=1
        if inputs['ggZZ']:  counter+=1
        if inputs['zjets']: counter+=1
        if inputs['ttbar']: counter+=1
        if inputs['zbb']:   counter+=1

        return counter

