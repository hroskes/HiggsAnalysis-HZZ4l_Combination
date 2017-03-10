#!/usr/bin/python
#-----------------------------------------------
# Latest update: 2012.08.30
# by Matt Snowball
#-----------------------------------------------
import sys, os, pwd, commands
import optparse, shlex, re
import math
from ROOT import *
import ROOT
from array import array
from properties_datacardClass import *
from properties_datacardClass_phase import *
from properties_datacardClass_2D import *
from inputReader import *

#define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-i', '--input', dest='inputDir', type='string', default="",    help='inputs directory')
    parser.add_option('-a', '--append', dest='appendName', type='string', default="",    help='append name for cards dir')
    parser.add_option('-b', action='store_true', dest='noX', default=True ,help='no X11 windows')
    parser.add_option('-t', '--templateDir', type='string', dest='templateDir', default="templates2D" ,help='directory with 2D template histos')
    parser.add_option('-m', '--model', type='string', dest='model', default="1D" ,help='model: 1D, phase, 2D ')
    parser.add_option('-r', '--datadir', type='string', dest='dataDirAppend', default="" ,help='dataDirAppend: Reference CMSdata folder per measurement')
    parser.add_option('-A', "--analysis", type='string', dest='analysis', default="fa3", help="fa3, fa2, or fL1")
    parser.add_option('-P', "--production", type='string', dest='production', help="see helperstuff/enums.py, search for class Production(MyEnum) for the list")
    parser.add_option('--mPOLE', '--mH', type='float', dest='mPOLE', default="125.6",
                      help='mPOLE (--mPOLE/--mH): Pole mass for the Higgs (def=125.6)')

    parser.add_option('--newMu', action='store_true', dest='newMu', default=False ,help='Use new muV=R*RV, muF=R*RF parameterization. Old approach is to have mu=r.')
    parser.add_option('--sigmaVVai', type='string', dest='sigmaVVai', default="",
                      help='sigmaVVai: Production (not *decay) sigma_ai(VBF+VH), follow T1,T2..T4..T7..T9, comma-separated (def=empty)')


    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()


    if (opt.appendName == ''):
        print 'Please pass an append name for the cards directory! Exiting...'
        sys.exit()

    if (opt.inputDir == ''):
        print 'Please pass an input directory! Exiting...'
        sys.exit()

    sigmaOpt = opt.sigmaVVai.split(",")
    opt.sigmaVVaiVal = dict()
    for i in range(1,10):
       opt.sigmaVVaiVal["T{0:.0f}".format(i)] = float(0.)
    if opt.newMu:
      for sopt in sigmaOpt:
        ssopt = sopt.split(":")
        if len(ssopt)==2:
          opt.sigmaVVaiVal[ssopt[0]]=float(ssopt[1])
          print "sigmaVVai for {0} = {1}".format(ssopt[0],opt.sigmaVVaiVal[ssopt[0]])



# define make directory function
def makeDirectory(subDirName):
    if (not os.path.exists(subDirName)):
        cmd = 'mkdir -p '+subDirName
        status, output = commands.getstatusoutput(cmd)
        if status !=0:
            print 'Error in creating submission dir '+subDirName+'. Exiting...'
            sys.exit()
    else:
        print 'Directory '+subDirName+' already exists. Exiting...'
    #    sys.exit()


#define function for processing of os command
def processCmd(cmd):
#    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()



def creationLoop(directory):
    global opt, args

    startMass=[ opt.mPOLE ]
    stepSizes=[ 0.1]
    endVal=[     1]

    if (opt.model== '1D'):
      myClass = properties_datacardClass(opt.analysis, opt.production)
    elif (opt.model== 'phase'):
      myClass = properties_datacardClass_phase()
    elif (opt.model== '2D'):
      myClass = properties_datacardClass_2D()
    else:
      print 'Unknown model '+opt.model+'. Exiting.'
      sys.exit()

    myClass.loadIncludes()

    myReader4e = inputReader(opt.inputDir+"/inputs_4e.txt")
    myReader4e.readInputs()
    theInputs4e = myReader4e.getInputs()

    myReader4mu = inputReader(opt.inputDir+"/inputs_4mu.txt")
    myReader4mu.readInputs()
    theInputs4mu = myReader4mu.getInputs()

    myReader2e2mu = inputReader(opt.inputDir+"/inputs_2e2mu.txt")
    myReader2e2mu.readInputs()
    theInputs2e2mu = myReader2e2mu.getInputs()


    a=0
    while (a < len(startMass) ):

	c = 0
        while (c < endVal[a] ):

            mStart = startMass[a]
            step = stepSizes[a]
            mh = mStart + ( step * c )
            mhs = str(mh).replace('.0','')

            print mh

            sqrts = str(theInputs4e['sqrts']).replace('.0','')

            makeDirectory(directory+'/HCG/'+mhs+'/'+sqrts+'TeV')
            makeDirectory(directory+'/HCG_XSxBR/'+mhs+'/'+sqrts+'TeV')
            myClass.makeCardsWorkspaces(mh,directory,theInputs4e,opt)
            myClass.makeCardsWorkspaces(mh,directory,theInputs4mu,opt)
            myClass.makeCardsWorkspaces(mh,directory,theInputs2e2mu,opt)

            c += 1


	a += 1






# the main procedure
def make_prop_DCsandWSs():

    # parse the arguments and options
    global opt, args
    parseOptions()

    if (opt.appendName != ''):
        dirName = 'cards_'+opt.appendName


    subdir = ['HCG','HCG_XSxBR','figs']

    for d in subdir:
        makeDirectory(dirName+'/'+d)


    creationLoop(dirName)


    sys.exit()



# run the create_RM_cfg() as main()
if __name__ == "__main__":
    make_prop_DCsandWSs()


