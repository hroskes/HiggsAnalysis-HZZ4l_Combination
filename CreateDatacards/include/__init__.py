import os
import ROOT
import subprocess

class cd(object):
    def __init__(self, folder):
        self.folder = folder
    def __enter__(self):
        self.bkpdir = os.getcwd()
        os.chdir(self.folder)
    def __exit__(self, *exceptioninfo):
        os.chdir(self.bkpdir)

includedir = os.path.relpath(os.path.dirname(__file__))

with cd(includedir):
    for f in os.listdir("."):
        if f.startswith("build"):
            subprocess.check_call(["root", "-l", "-b", "-q", f])


ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
ROOT.gSystem.AddIncludePath("-I"+os.path.join(includedir, '')) #-Iinclude/
ROOT.gROOT.ProcessLine(".L "+os.path.join(includedir, "tdrstyle.cc"))
ROOT.gSystem.Load("libRooFit")
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
ROOT.gSystem.Load(os.path.join(includedir, "HiggsCSandWidth_cc.so"))
ROOT.gSystem.Load(os.path.join(includedir, "HiggsCSandWidthSM4_cc.so"))
