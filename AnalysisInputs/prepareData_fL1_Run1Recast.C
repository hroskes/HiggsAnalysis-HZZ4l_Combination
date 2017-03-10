/* 
 * Prepare root files containing data events.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b prepareData.C+
 * This creates all files for 3 final states, 7 and 8 TeV and stores them in the final destination directory
 *
 */


//#define LINKMELA //Uncomment to link the MELA package to compute KD on the fly


#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TAxis.h"
#include "TFile.h"
#include "TLegend.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "TTree.h"
#include "TText.h"
#include "TStyle.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#ifdef LINKMELA
#include "ZZMatrixElement/MELA/interface/Mela.h"
#endif
#endif

using namespace std;

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//--- Flags to control re-computation of KD
bool recompute_ = false;       // Recompute KD (instead of taking it from the tree); when true, the following flags apply:
bool usePowhegTemplate=false;  // false use analytic bg
bool withPt_ = false;          // Include pT in KD
bool withY_  = false;          //    "    Y  "  "
int sqrts    = 8;              // sqrts, used only for withPt_/withY_

bool onlyICHEPStat = false;
bool dumpEvents = true;

#ifdef LINKMELA
Mela* myMELA;
#endif

int convertTreeForDatacards(TString inFile, TString outfile, bool useJET, bool VBFtag);

// Run all final states and sqrts in one go
void prepareData_fL1_Run1Recast() {

#ifdef LINKMELA
  if (recompute_) myMELA = new Mela(usePowhegTemplate); // this is safely leaked
#endif

  string appendDataName = "_fL1_Run1Recast_mH125";
  //  gSystem->Exec("mkdir -p "+ DataRootFilePath);
  gSystem->Exec("mkdir -p "+ DataRootFilePath + appendDataName);
  Int_t nOut[18];

  nOut[0]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu_Reprocessed_mH125.0.root", DataRootFilePath + appendDataName+"/hzz4mu_"  +lumistr7TeV+".root", false, false);
  nOut[1]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle_Reprocessed_mH125.0.root", DataRootFilePath + appendDataName+"/hzz4e_"   +lumistr7TeV+".root", false, false);
  nOut[2]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr_Reprocessed_mH125.0.root", DataRootFilePath + appendDataName+"/hzz2e2mu_"+lumistr7TeV+".root", false, false);
  nOut[3]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu_Reprocessed_mH125.0.root", DataRootFilePath + appendDataName+"/hzz4mu_"  +lumistr8TeV+".root", false, false);
  nOut[4]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle_Reprocessed_mH125.0.root", DataRootFilePath + appendDataName+"/hzz4e_"   +lumistr8TeV+".root", false, false);
  nOut[5]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr_Reprocessed_mH125.0.root", DataRootFilePath + appendDataName+"/hzz2e2mu_"+lumistr8TeV+".root", false, false);

  for (int ien = 0; ien<2; ien++){ for (int ich=0; ich<3; ich++) cout << "N Events: " << nOut[ien*3+ich] << endl; }
}


// The actual job
int convertTreeForDatacards(TString inFile, TString outfile, bool useJET, bool VBFtag){
  TChain* treedata;
  treedata= new TChain("SelectedTree");
  treedata->Add(inFile);

  int neventOut=0;
  Int_t run;
  Int_t ls=0;
  Long64_t event;
  int Lep1Id, Lep2Id, Lep3Id, Lep4Id;
  short Z1ids, Z2ids;

  int RunNumber=-1;
  Long64_t EventNumber=-1;

  float ZZMass;

  //Spin analysis variables
  float D_bkg_kin=-99;
  float D_g1q2=-99, D_g1q2int=-99;
  float D_g2=-99, D_g2int=-99, D_g2int_perp=-99;
  float D_g4=-99, D_g4int=-99, D_g4int_perp=-99;
  float D_ZG=-99, D_GG=-99;
  float D_ZG_PS=-99, D_GG_PS=-99;
  float D_ZGint=-99, D_GGint=-99;
  float D_ZG_PSint=-99, D_GG_PSint=-99;
  float D_ZG_L1=-99, D_ZG_L1int=-99, D_ZG_L1int_perp=-99;

  float D_bkg=-99;
  float D_bkg_prodIndep=-99;

  // Spin 1 and 2 discriminants
  float p1plusProdIndepKD = -99;
  float p1minusProdIndepKD = -99;
  float p1plusKD =-99;
  float p1minusKD =-99;
  float graviKD = -99;
  float qqgraviKD = -99;
  float p2h2plusKD = -99;
  float p2h2plus_qqb_KD = -99;
  float p2h3plusKD = -99;
  float p2h3plus_qqb_KD = -99;
  float p2hplusKD = -99;
  float p2hplus_qqb_KD = -99;
  float p2bplusKD = -99;
  float p2bplus_qqb_KD = -99;
  float p2h6plusKD = -99;
  float p2h6plus_qqb_KD = -99;
  float p2h7plusKD = -99;
  float p2h7plus_qqb_KD = -99;
  float p2hminusKD = -99;
  float p2hminus_qqb_KD = -99;
  float p2h9minusKD = -99;
  float p2h9minus_qqb_KD = -99;
  float p2h10minusKD = -99;
  float p2h10minus_qqb_KD = -99;
  float p2mProdIndepKD = -99;
  float p2h2plusProdIndepKD = -99;
  float p2h3plusProdIndepKD = -99;
  float p2hplusProdIndepKD = -99;
  float p2bplusProdIndepKD = -99;
  float p2h6plusProdIndepKD = -99;
  float p2h7plusProdIndepKD = -99;
  float p2hminusProdIndepKD = -99;
  float p2h9minusProdIndepKD = -99;
  float p2h10minusProdIndepKD = -99;

  treedata->SetBranchAddress("RunNumber", &RunNumber);
  treedata->SetBranchAddress("EventNumber", &EventNumber);

  treedata->SetBranchAddress("ZZMass", &ZZMass);

  treedata->SetBranchAddress("Z1ids", &Z1ids);
  treedata->SetBranchAddress("Z2ids", &Z2ids);
  treedata->SetBranchAddress("Lep1ID", &Lep1Id);
  treedata->SetBranchAddress("Lep2ID", &Lep2Id);
  treedata->SetBranchAddress("Lep3ID", &Lep3Id);
  treedata->SetBranchAddress("Lep4ID", &Lep4Id);

  treedata->SetBranchAddress("D_g1_vs_g2_phi0", &D_g2);
  treedata->SetBranchAddress("D_g1_vs_g4_phi0", &D_g4);
  treedata->SetBranchAddress("D_g2int_phi0", &D_g2int);
  treedata->SetBranchAddress("D_g4int_phi0", &D_g4int);
  treedata->SetBranchAddress("D_g2int_phi90", &D_g2int_perp);
  treedata->SetBranchAddress("D_g4int_phi90", &D_g4int_perp);
  treedata->SetBranchAddress("D_g1Q2_phi0", &D_g1q2);
  treedata->SetBranchAddress("D_g1Q2int_phi0", &D_g1q2int);

  treedata->SetBranchAddress("D_ZG", &D_ZG);
  treedata->SetBranchAddress("D_GG", &D_GG);
  treedata->SetBranchAddress("D_ZG_PS", &D_ZG_PS);
  treedata->SetBranchAddress("D_GG_PS", &D_GG_PS);
  treedata->SetBranchAddress("D_ZGint", &D_ZGint);
  treedata->SetBranchAddress("D_GGint", &D_GGint);
  treedata->SetBranchAddress("D_ZG_PSint", &D_ZG_PSint);
  treedata->SetBranchAddress("D_GG_PSint", &D_GG_PSint);
  treedata->SetBranchAddress("D_ZG_L1", &D_ZG_L1);
  treedata->SetBranchAddress("D_ZG_L1int_phi0", &D_ZG_L1int);
  treedata->SetBranchAddress("D_ZG_L1int_phi90", &D_ZG_L1int_perp);

  treedata->SetBranchAddress("D_bkg", &D_bkg);
  treedata->SetBranchAddress("D_bkg_kin", &D_bkg_kin);

  // Spin 1 and 2 discriminants
  treedata->SetBranchAddress("D_bkg_prodIndep", &D_bkg_prodIndep);

  treedata->SetBranchAddress("p1plusProdIndepKD", &p1plusProdIndepKD);
  treedata->SetBranchAddress("p1minusProdIndepKD", &p1minusProdIndepKD);
  treedata->SetBranchAddress("p1plusKD", &p1plusKD);
  treedata->SetBranchAddress("p1minusKD", &p1minusKD);

  treedata->SetBranchAddress("graviKD", &graviKD);
  treedata->SetBranchAddress("qqgraviKD", &qqgraviKD);
  treedata->SetBranchAddress("p2h2plusKD", &p2h2plusKD);
  treedata->SetBranchAddress("p2h2plus_qqb_KD", &p2h2plus_qqb_KD);
  treedata->SetBranchAddress("p2h3plusKD", &p2h3plusKD);
  treedata->SetBranchAddress("p2h3plus_qqb_KD", &p2h3plus_qqb_KD);
  treedata->SetBranchAddress("p2hplusKD", &p2hplusKD);
  treedata->SetBranchAddress("p2hplus_qqb_KD", &p2hplus_qqb_KD);
  treedata->SetBranchAddress("p2bplusKD", &p2bplusKD);
  treedata->SetBranchAddress("p2bplus_qqb_KD", &p2bplus_qqb_KD);
  treedata->SetBranchAddress("p2h6plusKD", &p2h6plusKD);
  treedata->SetBranchAddress("p2h6plus_qqb_KD", &p2h6plus_qqb_KD);
  treedata->SetBranchAddress("p2h7plusKD", &p2h7plusKD);
  treedata->SetBranchAddress("p2h7plus_qqb_KD", &p2h7plus_qqb_KD);
  treedata->SetBranchAddress("p2hminusKD", &p2hminusKD);
  treedata->SetBranchAddress("p2hminus_qqb_KD", &p2hminus_qqb_KD);
  treedata->SetBranchAddress("p2h9minusKD", &p2h9minusKD);
  treedata->SetBranchAddress("p2h9minus_qqb_KD", &p2h9minus_qqb_KD);
  treedata->SetBranchAddress("p2h10minusKD", &p2h10minusKD);
  treedata->SetBranchAddress("p2h10minus_qqb_KD", &p2h10minus_qqb_KD);
  treedata->SetBranchAddress("p2mProdIndepKD", &p2mProdIndepKD);
  treedata->SetBranchAddress("p2h2plusProdIndepKD", &p2h2plusProdIndepKD);
  treedata->SetBranchAddress("p2h3plusProdIndepKD", &p2h3plusProdIndepKD);
  treedata->SetBranchAddress("p2hplusProdIndepKD", &p2hplusProdIndepKD);
  treedata->SetBranchAddress("p2bplusProdIndepKD", &p2bplusProdIndepKD);
  treedata->SetBranchAddress("p2h6plusProdIndepKD", &p2h6plusProdIndepKD);
  treedata->SetBranchAddress("p2h7plusProdIndepKD", &p2h7plusProdIndepKD);
  treedata->SetBranchAddress("p2hminusProdIndepKD", &p2hminusProdIndepKD);
  treedata->SetBranchAddress("p2h9minusProdIndepKD", &p2h9minusProdIndepKD);
  treedata->SetBranchAddress("p2h10minusProdIndepKD", &p2h10minusProdIndepKD);

  TFile* newFile  = TFile::Open(outfile, "RECREATE");
  newFile->cd();
  TTree* newTree = new TTree("data_obs", "data_obs");
  Double_t CMS_zz4l_mass, CMS_zz4l_KD1, CMS_zz4l_KD2, CMS_zz4l_smd;

  newTree->Branch("CMS_zz4l_mass", &CMS_zz4l_mass, "CMS_zz4l_mass/D");
  //Spin Analysis Branches
  newTree->Branch("CMS_zz4l_smd", &CMS_zz4l_smd);
  newTree->Branch("CMS_zz4l_KD1", &CMS_zz4l_KD1);
  newTree->Branch("CMS_zz4l_KD2", &CMS_zz4l_KD2);

  newTree->Branch("CMS_zz4l_D_g1_vs_g2_phi0", &D_g2);
  newTree->Branch("CMS_zz4l_D_g1_vs_g4_phi0", &D_g4);
  newTree->Branch("CMS_zz4l_D_g2int_phi0", &D_g2int);
  newTree->Branch("CMS_zz4l_D_g4int_phi0", &D_g4int);
  newTree->Branch("CMS_zz4l_D_g2int_phi90", &D_g2int_perp);
  newTree->Branch("CMS_zz4l_D_g4int_phi90", &D_g4int_perp);
  newTree->Branch("CMS_zz4l_D_g1Q2_phi0", &D_g1q2);
  newTree->Branch("CMS_zz4l_D_g1Q2int_phi0", &D_g1q2int);

  newTree->Branch("CMS_zz4l_D_ZG", &D_ZG);
  newTree->Branch("CMS_zz4l_D_GG", &D_GG);
  newTree->Branch("CMS_zz4l_D_ZG_PS", &D_ZG_PS);
  newTree->Branch("CMS_zz4l_D_GG_PS", &D_GG_PS);
  newTree->Branch("CMS_zz4l_D_ZGint", &D_ZGint);
  newTree->Branch("CMS_zz4l_D_GGint", &D_GGint);
  newTree->Branch("CMS_zz4l_D_ZG_PSint", &D_ZG_PSint);
  newTree->Branch("CMS_zz4l_D_GG_PSint", &D_GG_PSint);
  newTree->Branch("CMS_zz4l_D_ZG_L1", &D_ZG_L1);
  newTree->Branch("CMS_zz4l_D_ZG_L1int_phi0", &D_ZG_L1int);
  newTree->Branch("CMS_zz4l_D_ZG_L1int_phi90", &D_ZG_L1int_perp);

  newTree->Branch("CMS_zz4l_D_bkg", &D_bkg);
  newTree->Branch("CMS_zz4l_D_bkg_kin", &D_bkg_kin);

  newTree->Branch("CMS_zz4l_p1plus_ProdIndepKD", &p1plusProdIndepKD);
  newTree->Branch("CMS_zz4l_p1minus_ProdIndepKD", &p1minusProdIndepKD);
  newTree->Branch("CMS_zz4l_p1plusKD", &p1plusKD);
  newTree->Branch("CMS_zz4l_p1minusKD", &p1minusKD);

  newTree->Branch("CMS_zz4l_graviKD", &graviKD);
  newTree->Branch("CMS_zz4l_qqgraviKD", &qqgraviKD);
  newTree->Branch("CMS_zz4l_p2h2plusKD", &p2h2plusKD);
  newTree->Branch("CMS_zz4l_p2h2plus_qqb_KD", &p2h2plus_qqb_KD);
  newTree->Branch("CMS_zz4l_p2h3plusKD", &p2h3plusKD);
  newTree->Branch("CMS_zz4l_p2h3plus_qqb_KD", &p2h3plus_qqb_KD);
  newTree->Branch("CMS_zz4l_p2hplusKD", &p2hplusKD);
  newTree->Branch("CMS_zz4l_p2hplus_qqb_KD", &p2hplus_qqb_KD);
  newTree->Branch("CMS_zz4l_p2bplusKD", &p2bplusKD);
  newTree->Branch("CMS_zz4l_p2bplus_qqb_KD", &p2bplus_qqb_KD);
  newTree->Branch("CMS_zz4l_p2h6plusKD", &p2h6plusKD);
  newTree->Branch("CMS_zz4l_p2h6plus_qqb_KD", &p2h6plus_qqb_KD);
  newTree->Branch("CMS_zz4l_p2h7plusKD", &p2h7plusKD);
  newTree->Branch("CMS_zz4l_p2h7plus_qqb_KD", &p2h7plus_qqb_KD);
  newTree->Branch("CMS_zz4l_p2hminusKD", &p2hminusKD);
  newTree->Branch("CMS_zz4l_p2hminus_qqb_KD", &p2hminus_qqb_KD);
  newTree->Branch("CMS_zz4l_p2h9minusKD", &p2h9minusKD);
  newTree->Branch("CMS_zz4l_p2h9minus_qqb_KD", &p2h9minus_qqb_KD);
  newTree->Branch("CMS_zz4l_p2h10minusKD", &p2h10minusKD);
  newTree->Branch("CMS_zz4l_p2h10minus_qqb_KD", &p2h10minus_qqb_KD);
  newTree->Branch("CMS_zz4l_p2mProdIndepKD", &p2mProdIndepKD);
  newTree->Branch("CMS_zz4l_p2h2plusProdIndepKD", &p2h2plusProdIndepKD);
  newTree->Branch("CMS_zz4l_p2h3plusProdIndepKD", &p2h3plusProdIndepKD);
  newTree->Branch("CMS_zz4l_p2hplusProdIndepKD", &p2hplusProdIndepKD);
  newTree->Branch("CMS_zz4l_p2bplusProdIndepKD", &p2bplusProdIndepKD);
  newTree->Branch("CMS_zz4l_p2h6plusProdIndepKD", &p2h6plusProdIndepKD);
  newTree->Branch("CMS_zz4l_p2h7plusProdIndepKD", &p2h7plusProdIndepKD);
  newTree->Branch("CMS_zz4l_p2hminusProdIndepKD", &p2hminusProdIndepKD);
  newTree->Branch("CMS_zz4l_p2h9minusProdIndepKD", &p2h9minusProdIndepKD);
  newTree->Branch("CMS_zz4l_p2h10minusProdIndepKD", &p2h10minusProdIndepKD);


  cout << inFile << " entries: " << treedata->GetEntries() << endl;
  for (int iEvt=0; iEvt<treedata->GetEntries(); iEvt++){
    treedata->GetEntry(iEvt);

    if (Z1ids==-13*13){
      Lep1Id-13;
      Lep2Id=-13;
    }
    else if (Z1ids==-11*11){
      Lep1Id=11;
      Lep2Id=-11;
    }
    if (Z2ids==-13*13){
      Lep3Id=13;
      Lep4Id=-13;
    }
    else if (Z2ids==-11*11){
      Lep3Id=11;
      Lep4Id=-11;
    }
    int lepIdOrdered[4]={ Lep1Id, Lep2Id, Lep3Id, Lep4Id };

    CMS_zz4l_mass = double(ZZMass);
    CMS_zz4l_smd = D_bkg;
    CMS_zz4l_KD1 = D_g1q2;
    CMS_zz4l_KD2 = D_g2;

    ++neventOut;
    newTree->Fill();
    if (dumpEvents) {
      if (!useJET && !VBFtag) {
        cout.setf(ios::fixed);
        int ifs=-1;
        if (outfile.Contains("4e")) ifs = 0;
        else if (outfile.Contains("4mu")) ifs = 1;
        else if (outfile.Contains("2e2mu")) ifs = 2;
        cout << run << ":" << ls << ":" << event << " " << ifs << " " << CMS_zz4l_mass << " " << CMS_zz4l_smd << " " << CMS_zz4l_KD1 << " " << CMS_zz4l_KD2 << endl;
      }
    }
  }
  newTree->Write("data_obs");
  newFile->Close();

  cout << "written: " << outfile << " entries: " << neventOut << endl << endl;
  return  neventOut;

}

