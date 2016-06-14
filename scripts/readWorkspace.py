#!/usr/bin/env python
import os
import ROOT
import sys

ROOT.gROOT.LoadMacro(os.path.join(os.environ["CMSSW_BASE"], "src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards/readWorkspace.c"))
ROOT.readWorkspace(sys.argv[1])
