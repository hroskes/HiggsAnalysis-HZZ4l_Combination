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
void prepareData() {

#ifdef LINKMELA
  if (recompute_) myMELA = new Mela(usePowhegTemplate); // this is safely leaked
#endif

  gSystem->Exec("mkdir -p "+ DataRootFilePath);
//  gSystem->Exec("mkdir -p "+ DataRootFilePath + "_fL1");
//  gSystem->Exec("mkdir -p "+ DataRootFilePath + "_fZgsL1");
//  gSystem->Exec("mkdir -p "+ DataRootFilePath + "_fa3perp");
  Int_t nOut[18];

  nOut[0]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"/hzz4mu_"  +lumistr7TeV+"_1.root",true, true);
  nOut[1]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"/hzz4e_"   +lumistr7TeV+"_1.root",true, true);
  nOut[2]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"/hzz2e2mu_"+lumistr7TeV+"_1.root",true, true);
  nOut[3]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"/hzz4mu_"  +lumistr8TeV+"_1.root",true, true);
  nOut[4]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"/hzz4e_"   +lumistr8TeV+"_1.root",true, true);
  nOut[5]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"/hzz2e2mu_"+lumistr8TeV+"_1.root",true, true);

  nOut[6]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"/hzz4mu_"  +lumistr7TeV+"_0.root",true, false);
  nOut[7]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"/hzz4e_"   +lumistr7TeV+"_0.root",true, false);
  nOut[8]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"/hzz2e2mu_"+lumistr7TeV+"_0.root",true, false);
  nOut[9]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"/hzz4mu_"  +lumistr8TeV+"_0.root",true, false);
  nOut[10]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"/hzz4e_"   +lumistr8TeV+"_0.root",true, false);
  nOut[11]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"/hzz2e2mu_"+lumistr8TeV+"_0.root",true, false);

  nOut[12]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"/hzz4mu_"  +lumistr7TeV+".root",false, false);
  nOut[13]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"/hzz4e_"   +lumistr7TeV+".root",false, false);
  nOut[14]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"/hzz2e2mu_"+lumistr7TeV+".root",false, false);
  nOut[15]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"/hzz4mu_"  +lumistr8TeV+".root",false, false);
  nOut[16]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"/hzz4e_"   +lumistr8TeV+".root",false, false);
  nOut[17]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"/hzz2e2mu_"+lumistr8TeV+".root",false, false);

  TString jetsT[3]={"DiJet     ",
		    "Untagged  ",
		    "Inclusive "};
  TString energies[2]={"7TeV", "8TeV "};
  TString channels[3]={"4mu   ",
		       "4e    ",
		       "2e2mu "};
  cout<<"Summary  "<<endl;
  for(int ijet=0;ijet<3;ijet++)
    for(int ien = 0;ien<2;ien++)
      for(int ich=0;ich<3;ich++)
	cout<<"N Events out for "<<jetsT[ijet]<<energies[ien].Data()<<channels[ich].Data()<<": "<<nOut[ijet*6+ien*3+ich]<<endl;
}


// The actual job
int convertTreeForDatacards(TString inFile, TString outfile, bool useJET, bool VBFtag){
  TChain* treedata ;
  treedata= new TChain("SelectedTree");
  treedata->Add(inFile);

  int neventOut=0;
  Int_t run;
  Int_t ls=0;
  Long64_t event;
		int Lep1Id,Lep2Id,Lep3Id,Lep4Id;
		short Z1ids,Z2ids;
  float mzz, mzzErr;
  float m1, m2;
  //float ZZVAKD;
  float pt4l, fisher;
  int NJets_30;

  //Spin analysis variables
		float p0plus_VAJHU;
		float p0hplus_VAJHU;
		float p0minus_VAJHU;
		float p0_g1prime2_VAJHU;
		float pg1g1prime2_VAJHU;
		float pg1g2_VAJHU;
		float pg1g4_VAJHU;
		float pg1g2_pi2_VAJHU;
		float pg1g4_pi2_VAJHU;

		float pzzzg_VAJHU;
		float pzzgg_VAJHU;
		float pzzzg_PS_VAJHU;
		float pzzgg_PS_VAJHU;
		float p0Zgs_VAJHU;
		float p0gsgs_VAJHU;
		float p0Zgs_PS_VAJHU;
		float p0gsgs_PS_VAJHU;
		float p0Zgs_g1prime2_VAJHU;
		float pzzzgs_g1prime2_VAJHU;
		float pzzzgs_g1prime2_pi2_VAJHU;
		
		float p1plusProdIndepVA, p1minusProdIndepVA;
		float p1plus_VAJHU,p1_VAJHU;
		float p2minimalVA, p2minimalVA_qqb, p2h2plusVA, p2h2plusVA_qqb;
		float p2h3plusVA, p2h3plusVA_qqb, p2hplusVA, p2hplusVA_qqb;
		float p2bplusVA, p2bplusVA_qqb, p2h6plusVA, p2h6plusVA_qqb, p2h7plusVA, p2h7plusVA_qqb, p2hminusVA, p2hminusVA_qqb, p2h9minusVA, p2h9minusVA_qqb, p2h10minusVA, p2h10minusVA_qqb;
		float p2mProdIndepVA, p2h2plusProdIndepVA, p2h3plusProdIndepVA, p2hplusProdIndepVA, p2bplusProdIndepVA, p2h6plusProdIndepVA, p2h7plusProdIndepVA, p2hminusProdIndepVA, p2h9minusProdIndepVA, p2h10minusProdIndepVA;

		float bkg_VAMCFM, bkg_prodIndep_VAMCFM;
		float bkg_VAMCFM_OLD,ggzz_p0plus_VAMCFM_OLD,ggzz_VAMCFM_OLD,p0plus_VAMCFM_OLD;
		float ggzz_VAMCFM,ggzz_c1_VAMCFM,ggzz_c5_VAMCFM;
		float bkg_VAMCFM_noscale , ggzz_VAMCFM_noscale , ggHZZ_prob_pure , ggHZZ_prob_int , ggHZZ_prob_pure_noscale , ggHZZ_prob_int_noscale;
		float qqScale , ggScale;
		float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
		float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;


		float Dgg10_VAMCFM=-99;

		double D_ggZZ=-99;
		double D_bkg_kin=-99;
		double D_Gamma_r1=-99;
		double D_Gamma_r5=-99;
		double D_Gamma_r10=-99;
		double D_Gamma_r15=-99;
		double D_Gamma_r20=-99;
		double D_Gamma_r25=-99;
		double D_Gamma_gg_r1=-99;
		double D_Gamma_gg_r5=-99;
		double D_Gamma_gg_r10=-99;
		double D_Gamma_gg_r15=-99;
		double D_Gamma_gg_r20=-99;
		double D_Gamma_gg_r25=-99;
		double D_Gamma=-99;
		double D_Gamma_int=-99;
		double D_g1q2=-99,D_g1q2int=-99;
		double D_g2=-99,D_g2int=-99,D_g2int_perp=-99;
		double D_g4=-99,D_g4int=-99,D_g4int_perp=-99;
		double D_ZG=-99, D_GG=-99;
		double D_ZG_PS=-99, D_GG_PS=-99;
		double D_ZGint=-99, D_GGint=-99;
		double D_ZG_PSint=-99, D_GG_PSint=-99;
		double D_ZG_L1=-99,D_ZG_L1int=-99,D_ZG_L1int_perp=-99;

		double D_bkg=-99;
		double D_bkg_prodIndep=-99;
		double D_bkg_ScaleUp=-99;
		double D_bkg_prodIndep_ScaleUp=-99;
		double D_bkg_ScaleDown=-99;
		double D_bkg_prodIndep_ScaleDown=-99;
		double D_bkg_ResUp=-99;
		double D_bkg_prodIndep_ResUp=-99;
		double D_bkg_ResDown=-99;
		double D_bkg_prodIndep_ResDown=-99;

// Spin 1 and 2 discriminants
		double p1plusProdIndepKD = -99;
		double p1minusProdIndepKD = -99;
		double p1plusKD =-99;
		double p1minusKD =-99;
		double graviKD = -99;
		double qqgraviKD = -99;
		double p2h2plusKD = -99;
		double p2h2plus_qqb_KD = -99;
		double p2h3plusKD = -99;
		double p2h3plus_qqb_KD = -99;
		double p2hplusKD = -99;
		double p2hplus_qqb_KD = -99;
		double p2bplusKD = -99;
		double p2bplus_qqb_KD = -99;
		double p2h6plusKD = -99;
		double p2h6plus_qqb_KD = -99;
		double p2h7plusKD = -99;
		double p2h7plus_qqb_KD = -99;
		double p2hminusKD = -99;
		double p2hminus_qqb_KD = -99;
		double p2h9minusKD = -99;
		double p2h9minus_qqb_KD = -99;
		double p2h10minusKD = -99;
		double p2h10minus_qqb_KD = -99;
		double p2mProdIndepKD = -99;
		double p2h2plusProdIndepKD = -99;
		double p2h3plusProdIndepKD = -99;
		double p2hplusProdIndepKD = -99;
		double p2bplusProdIndepKD = -99;
		double p2h6plusProdIndepKD = -99;
		double p2h7plusProdIndepKD = -99;
		double p2hminusProdIndepKD = -99;
		double p2h9minusProdIndepKD = -99;
		double p2h10minusProdIndepKD = -99;

  treedata->SetBranchAddress("RunNumber",&run);
  //  treedata->SetBranchAddress("LumiNumber",&ls);
  treedata->SetBranchAddress("EventNumber",&event);
  treedata->SetBranchAddress("ZZMass",&mzz);
  treedata->SetBranchAddress("ZZMassErrCorr",&mzzErr);
  treedata->SetBranchAddress("Z1Mass",&m1);
  treedata->SetBranchAddress("Z2Mass",&m2);
			treedata->SetBranchAddress("Z1ids", &Z1ids);
			treedata->SetBranchAddress("Z2ids", &Z2ids);
			treedata->SetBranchAddress("Lep1ID", &Lep1Id);
			treedata->SetBranchAddress("Lep2ID", &Lep2Id);
			treedata->SetBranchAddress("Lep3ID", &Lep3Id);
			treedata->SetBranchAddress("Lep4ID", &Lep4Id);

  treedata->SetBranchAddress("ZZPt",&pt4l);
  //treedata->SetBranchAddress("ZZRapidity",&Y4l);
  treedata->SetBranchAddress("ZZFisher",&fisher);

  treedata->SetBranchAddress("NJets30", &NJets_30);
  treedata->SetBranchAddress("Dgg10_VAMCFM",&Dgg10_VAMCFM);

  //Spin Disc variables
			treedata->SetBranchAddress("p0plus_VAJHU",&p0plus_VAJHU);
			treedata->SetBranchAddress("p0hplus_VAJHU",&p0hplus_VAJHU);
			treedata->SetBranchAddress("p0minus_VAJHU",&p0minus_VAJHU);
			treedata->SetBranchAddress("p0_g1prime2_VAJHU",&p0_g1prime2_VAJHU);
			treedata->SetBranchAddress("pg1g1prime2_VAJHU",&pg1g1prime2_VAJHU);
			treedata->SetBranchAddress("pg1g2_VAJHU",&pg1g2_VAJHU);
			treedata->SetBranchAddress("pg1g4_VAJHU",&pg1g4_VAJHU);
			treedata->SetBranchAddress("pg1g2_pi2_VAJHU",&pg1g2_pi2_VAJHU);
			treedata->SetBranchAddress("pg1g4_pi2_VAJHU",&pg1g4_pi2_VAJHU);

			treedata->SetBranchAddress("pzzzg_VAJHU",&pzzzg_VAJHU);
			treedata->SetBranchAddress("pzzgg_VAJHU",&pzzgg_VAJHU);
			treedata->SetBranchAddress("pzzzg_PS_VAJHU",&pzzzg_PS_VAJHU);
			treedata->SetBranchAddress("pzzgg_PS_VAJHU",&pzzgg_PS_VAJHU);
			treedata->SetBranchAddress("p0Zgs_VAJHU",&p0Zgs_VAJHU);
			treedata->SetBranchAddress("p0gsgs_VAJHU",&p0gsgs_VAJHU);
			treedata->SetBranchAddress("p0Zgs_PS_VAJHU",&p0Zgs_PS_VAJHU);
			treedata->SetBranchAddress("p0gsgs_PS_VAJHU",&p0gsgs_PS_VAJHU);
			treedata->SetBranchAddress("p0Zgs_g1prime2_VAJHU",&p0Zgs_g1prime2_VAJHU);
			treedata->SetBranchAddress("pzzzgs_g1prime2_VAJHU",&pzzzgs_g1prime2_VAJHU);
			treedata->SetBranchAddress("pzzzgs_g1prime2_pi2_VAJHU",&pzzzgs_g1prime2_pi2_VAJHU);

			treedata->SetBranchAddress("p0plus_VAMCFM",&p0plus_VAMCFM_OLD);
			treedata->SetBranchAddress("bkg_VAMCFM",&bkg_VAMCFM_OLD);
			treedata->SetBranchAddress("ggzz_p0plus_VAMCFM",&ggzz_p0plus_VAMCFM_OLD);
			treedata->SetBranchAddress("ggzz_VAMCFM",&ggzz_VAMCFM_OLD);

			treedata->SetBranchAddress("p0plus_m4l",&p0plus_m4l);
			treedata->SetBranchAddress("bkg_m4l",&bkg_m4l);
			treedata->SetBranchAddress("p0plus_m4l_ScaleUp",&p0plus_m4l_ScaleUp);
			treedata->SetBranchAddress("bkg_m4l_ScaleUp",&bkg_m4l_ScaleUp);
			treedata->SetBranchAddress("p0plus_m4l_ScaleDown",&p0plus_m4l_ScaleDown);
			treedata->SetBranchAddress("bkg_m4l_ScaleDown",&bkg_m4l_ScaleDown);
			treedata->SetBranchAddress("p0plus_m4l_ResUp",&p0plus_m4l_ResUp);
			treedata->SetBranchAddress("bkg_m4l_ResUp",&bkg_m4l_ResUp);
			treedata->SetBranchAddress("p0plus_m4l_ResDown",&p0plus_m4l_ResDown);
			treedata->SetBranchAddress("bkg_m4l_ResDown",&bkg_m4l_ResDown);

			treedata->SetBranchAddress("p1plus_prodIndep_VAJHU", &p1plusProdIndepVA);
			treedata->SetBranchAddress("p1_prodIndep_VAJHU", &p1minusProdIndepVA);
			treedata->SetBranchAddress("p1plus_VAJHU", &p1plus_VAJHU);
			treedata->SetBranchAddress("p1_VAJHU", &p1_VAJHU);

			treedata->SetBranchAddress("p2_VAJHU", &p2minimalVA);
			treedata->SetBranchAddress("p2qqb_VAJHU", &p2minimalVA_qqb);
			treedata->SetBranchAddress("p2h2plus_gg_VAJHU", &p2h2plusVA);
			treedata->SetBranchAddress("p2h2plus_qqbar_VAJHU", &p2h2plusVA_qqb);
			treedata->SetBranchAddress("p2h3plus_gg_VAJHU", &p2h3plusVA);
			treedata->SetBranchAddress("p2h3plus_qqbar_VAJHU", &p2h3plusVA_qqb);
			treedata->SetBranchAddress("p2hplus_VAJHU", &p2hplusVA);
			treedata->SetBranchAddress("p2hplus_qqb_VAJHU", &p2hplusVA_qqb);
			treedata->SetBranchAddress("p2bplus_VAJHU", &p2bplusVA);
			treedata->SetBranchAddress("p2bplus_qqb_VAJHU", &p2bplusVA_qqb);
			treedata->SetBranchAddress("p2h6plus_gg_VAJHU", &p2h6plusVA);
			treedata->SetBranchAddress("p2h6plus_qqbar_VAJHU", &p2h6plusVA_qqb);
			treedata->SetBranchAddress("p2h7plus_gg_VAJHU", &p2h7plusVA);
			treedata->SetBranchAddress("p2h7plus_qqbar_VAJHU", &p2h7plusVA_qqb);
			treedata->SetBranchAddress("p2hminus_VAJHU", &p2hminusVA);
			treedata->SetBranchAddress("p2hminus_qqb_VAJHU", &p2hminusVA_qqb);
			treedata->SetBranchAddress("p2h9minus_gg_VAJHU", &p2h9minusVA);
			treedata->SetBranchAddress("p2h9minus_qqbar_VAJHU", &p2h9minusVA_qqb);
			treedata->SetBranchAddress("p2h10minus_gg_VAJHU", &p2h10minusVA);
			treedata->SetBranchAddress("p2h10minus_qqbar_VAJHU", &p2h10minusVA_qqb);
			treedata->SetBranchAddress("p2_prodIndep_VAJHU", &p2mProdIndepVA);
			treedata->SetBranchAddress("p2h2plus_prodIndep_VAJHU", &p2h2plusProdIndepVA);
			treedata->SetBranchAddress("p2h3plus_prodIndep_VAJHU", &p2h3plusProdIndepVA);
			treedata->SetBranchAddress("p2hplus_prodIndep_VAJHU", &p2hplusProdIndepVA);
			treedata->SetBranchAddress("p2bplus_prodIndep_VAJHU", &p2bplusProdIndepVA);
			treedata->SetBranchAddress("p2h6plus_prodIndep_VAJHU", &p2h6plusProdIndepVA);
			treedata->SetBranchAddress("p2h7plus_prodIndep_VAJHU", &p2h7plusProdIndepVA);
			treedata->SetBranchAddress("p2hminus_prodIndep_VAJHU", &p2hminusProdIndepVA);
			treedata->SetBranchAddress("p2h9minus_prodIndep_VAJHU", &p2h9minusProdIndepVA);
			treedata->SetBranchAddress("p2h10minus_prodIndep_VAJHU", &p2h10minusProdIndepVA);
			treedata->SetBranchAddress("bkg_prodIndep_VAMCFM",&bkg_prodIndep_VAMCFM);
 
   
  TFile* newFile  = new TFile(outfile, "RECREATE");
  newFile->cd();
  TTree* newTree = new TTree("data_obs","data_obs"); 
  Double_t CMS_zz4l_mass, melaLD, CMS_zz4l_massErr, CMS_zz4l_massRelErr;
  Double_t pt = -99, Fisher = -99;

  newTree->Branch("CMS_zz4l_mass",&CMS_zz4l_mass,"CMS_zz4l_mass/D");
  newTree->Branch("CMS_zz4l_massErr",&CMS_zz4l_massErr,"CMS_zz4l_massErr/D");
  newTree->Branch("CMS_zz4l_massRelErr",&CMS_zz4l_massRelErr,"CMS_zz4l_massRelErr/D");
  newTree->Branch("melaLD",&melaLD,"melaLD/D");
  newTree->Branch("CMS_zz4l_Fisher",&Fisher,"CMS_zz4l_Fisher/D");
  newTree->Branch("CMS_zz4l_Pt",&pt,"CMS_zz4l_Pt/D");

  //Spin Analysis Branches
//			newTree->Branch("CMS_zz4l_D_g1_vs_g2_phi0",&D_g2);
			newTree->Branch("CMS_zz4l_KD2",&D_g2);
			newTree->Branch("CMS_zz4l_D_g1_vs_g4_phi0",&D_g4);
			newTree->Branch("CMS_zz4l_D_g2int_phi0",&D_g2int);
			newTree->Branch("CMS_zz4l_D_g4int_phi0",&D_g4int);
			newTree->Branch("CMS_zz4l_D_g2int_phi90",&D_g2int_perp);
			newTree->Branch("CMS_zz4l_D_g4int_phi90",&D_g4int_perp);
			newTree->Branch("CMS_zz4l_KD1",&D_g1q2);
//			newTree->Branch("CMS_zz4l_D_g1Q2_phi0",&D_g1q2);
			newTree->Branch("CMS_zz4l_D_g1Q2int_phi0",&D_g1q2int);

			newTree->Branch("CMS_zz4l_D_ZG",&D_ZG);
			newTree->Branch("CMS_zz4l_D_GG",&D_GG);
			newTree->Branch("CMS_zz4l_D_ZG_PS",&D_ZG_PS);
			newTree->Branch("CMS_zz4l_D_GG_PS",&D_GG_PS);
			newTree->Branch("CMS_zz4l_D_ZGint",&D_ZGint);
			newTree->Branch("CMS_zz4l_D_GGint",&D_GGint);
			newTree->Branch("CMS_zz4l_D_ZG_PSint",&D_ZG_PSint);
			newTree->Branch("CMS_zz4l_D_GG_PSint",&D_GG_PSint);
			newTree->Branch("D_ZG_L1",&D_ZG_L1);
			newTree->Branch("D_ZG_L1int_phi0",&D_ZG_L1int);
			newTree->Branch("D_ZG_L1int_phi90",&D_ZG_L1int_perp);

			newTree->Branch("CMS_zz4l_smd",&D_bkg);
			newTree->Branch("CMS_zz4l_D_bkg_kin",&D_bkg_kin);
			newTree->Branch("CMS_zz4l_D_bkg_ScaleUp",&D_bkg_ScaleUp);
			newTree->Branch("CMS_zz4l_D_bkg_ScaleDown",&D_bkg_ScaleDown);
			newTree->Branch("CMS_zz4l_D_bkg_ResUp",&D_bkg_ResUp);
			newTree->Branch("CMS_zz4l_D_bkg_ResDown",&D_bkg_ResDown);

			newTree->Branch("CMS_zz4l_Dgg10_VAMCFM",&D_Gamma_r10);

			newTree->Branch("CMS_zz4l_D_bkg_prodIndep",&D_bkg_prodIndep);
			newTree->Branch("CMS_zz4l_D_bkg_prodIndep_ScaleUp",&D_bkg_prodIndep_ScaleUp);
			newTree->Branch("CMS_zz4l_D_bkg_prodIndep_ScaleDown",&D_bkg_prodIndep_ScaleDown);
			newTree->Branch("CMS_zz4l_D_bkg_prodIndep_ResUp",&D_bkg_prodIndep_ResUp);
			newTree->Branch("CMS_zz4l_D_bkg_prodIndep_ResDown",&D_bkg_prodIndep_ResDown);

			newTree->Branch("CMS_zz4l_p1plus_ProdIndepKD", &p1plusProdIndepKD);
			newTree->Branch("CMS_zz4l_p1minus_ProdIndepKD", &p1minusProdIndepKD);
			newTree->Branch("CMS_zz4l_p1plusKD", &p1plusKD);
			newTree->Branch("CMS_zz4l_p1minusKD", &p1minusKD);

			newTree->Branch("CMS_zz4l_graviKD",&graviKD);
			newTree->Branch("CMS_zz4l_qqgraviKD",&qqgraviKD);
			newTree->Branch("CMS_zz4l_p2h2plusKD",&p2h2plusKD);
			newTree->Branch("CMS_zz4l_p2h2plus_qqb_KD",&p2h2plus_qqb_KD);
			newTree->Branch("CMS_zz4l_p2h3plusKD",&p2h3plusKD);
			newTree->Branch("CMS_zz4l_p2h3plus_qqb_KD",&p2h3plus_qqb_KD);
			newTree->Branch("CMS_zz4l_p2hplusKD",&p2hplusKD);
			newTree->Branch("CMS_zz4l_p2hplus_qqb_KD",&p2hplus_qqb_KD);
			newTree->Branch("CMS_zz4l_p2bplusKD",&p2bplusKD);
			newTree->Branch("CMS_zz4l_p2bplus_qqb_KD",&p2bplus_qqb_KD);
			newTree->Branch("CMS_zz4l_p2h6plusKD",&p2h6plusKD);
			newTree->Branch("CMS_zz4l_p2h6plus_qqb_KD",&p2h6plus_qqb_KD);
			newTree->Branch("CMS_zz4l_p2h7plusKD",&p2h7plusKD);
			newTree->Branch("CMS_zz4l_p2h7plus_qqb_KD",&p2h7plus_qqb_KD);
			newTree->Branch("CMS_zz4l_p2hminusKD",&p2hminusKD);
			newTree->Branch("CMS_zz4l_p2hminus_qqb_KD",&p2hminus_qqb_KD);
			newTree->Branch("CMS_zz4l_p2h9minusKD",&p2h9minusKD);
			newTree->Branch("CMS_zz4l_p2h9minus_qqb_KD",&p2h9minus_qqb_KD);
			newTree->Branch("CMS_zz4l_p2h10minusKD",&p2h10minusKD);
			newTree->Branch("CMS_zz4l_p2h10minus_qqb_KD",&p2h10minus_qqb_KD);
			newTree->Branch("CMS_zz4l_p2mProdIndepKD",&p2mProdIndepKD);
			newTree->Branch("CMS_zz4l_p2h2plusProdIndepKD",&p2h2plusProdIndepKD);
			newTree->Branch("CMS_zz4l_p2h3plusProdIndepKD",&p2h3plusProdIndepKD);
			newTree->Branch("CMS_zz4l_p2hplusProdIndepKD",&p2hplusProdIndepKD);
			newTree->Branch("CMS_zz4l_p2bplusProdIndepKD",&p2bplusProdIndepKD);
			newTree->Branch("CMS_zz4l_p2h6plusProdIndepKD",&p2h6plusProdIndepKD);
			newTree->Branch("CMS_zz4l_p2h7plusProdIndepKD",&p2h7plusProdIndepKD);
			newTree->Branch("CMS_zz4l_p2hminusProdIndepKD",&p2hminusProdIndepKD);
			newTree->Branch("CMS_zz4l_p2h9minusProdIndepKD",&p2h9minusProdIndepKD);
			newTree->Branch("CMS_zz4l_p2h10minusProdIndepKD",&p2h10minusProdIndepKD);


  cout << inFile << " entries: " << treedata->GetEntries() << endl;
  for(int iEvt=0; iEvt<treedata->GetEntries(); iEvt++){
    //    if(iEvt%5000==0) cout << "event: " << iEvt << endl;
    treedata->GetEntry(iEvt);

    //cout << run << endl;
    
    if (onlyICHEPStat && run>=198049) continue;

    //if( NJets_30 > 0) cout << "NJets_30: " << NJets_30 << endl;

    if ((useJET && VBFtag && NJets_30 < 2) || (useJET && !VBFtag && NJets_30 >= 2)) continue;

				if(Z1ids==-13*13){
					Lep1Id=-13;
					Lep2Id=13;
				}
				else if(Z1ids==-11*11){
					Lep1Id=-11;
					Lep2Id=11;
				};
				if(Z2ids==-13*13){
					Lep3Id=-13;
					Lep4Id=13;
				}
				else if(Z2ids==-11*11){
					Lep3Id=-11;
					Lep4Id=11;
				};
			int lepIdOrdered[4]={ Lep1Id,Lep2Id,Lep3Id,Lep4Id };

    CMS_zz4l_mass = double(mzz);
    CMS_zz4l_massErr = double(mzzErr);
    CMS_zz4l_massRelErr = double(mzzErr/mzz);
    //pseudomelaLD = pseudomela;
    //melaLD = mela;
    //supermelaLD = 0;
    melaLD=double(p0plus_VAJHU/(p0plus_VAJHU+bkg_VAMCFM));
    //ZZVAKD=double(p0plus_VAJHU/(bkg_VAMCFM + p0plus_VAJHU));
    if(useJET && !VBFtag) 
      {
	if(pt4l > 200.)
	  {
	    //if pt out of range set to middle of higest bin
	    pt4l = 198.;
	  }
	pt = double(pt4l);
      }
    if(useJET && VBFtag) 
      {
	if(fisher > 2.)
	  {
	    //if fisher out of range set to middle of higest bin
	    fisher = 1.98;
	  }
	Fisher = double(fisher);
      }

    //make variables for spin analysis

			bkg_VAMCFM=bkg_VAMCFM_OLD;
			ggzz_VAMCFM_noscale=ggzz_VAMCFM_OLD;
			ggHZZ_prob_pure_noscale = p0plus_VAMCFM_OLD;
			ggHZZ_prob_int_noscale = ggzz_p0plus_VAMCFM_OLD;

			ggHZZ_prob_int_noscale = ggHZZ_prob_int_noscale - ggHZZ_prob_pure_noscale - ggzz_VAMCFM_noscale;

			D_bkg_kin = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM );
			D_bkg = p0plus_VAJHU*p0plus_m4l / ( p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l );
			D_bkg_ScaleUp = (p0plus_VAJHU*p0plus_m4l_ScaleUp) / (p0plus_VAJHU*p0plus_m4l_ScaleUp + bkg_VAMCFM_OLD*bkg_m4l_ScaleUp);
			D_bkg_ScaleDown = (p0plus_VAJHU*p0plus_m4l_ScaleDown) / (p0plus_VAJHU*p0plus_m4l_ScaleDown + bkg_VAMCFM_OLD*bkg_m4l_ScaleDown);
			D_bkg_ResUp = (p0plus_VAJHU*p0plus_m4l_ResUp) / (p0plus_VAJHU*p0plus_m4l_ResUp + bkg_VAMCFM_OLD*bkg_m4l_ResUp);
			D_bkg_ResDown = (p0plus_VAJHU*p0plus_m4l_ResDown) / (p0plus_VAJHU*p0plus_m4l_ResDown + bkg_VAMCFM_OLD*bkg_m4l_ResDown);

				D_bkg_prodIndep = (p0plus_VAJHU*p0plus_m4l)/(p0plus_VAJHU*p0plus_m4l + bkg_prodIndep_VAMCFM*bkg_m4l);
				D_bkg_prodIndep_ScaleUp   = (p0plus_VAJHU*p0plus_m4l_ScaleUp) / (p0plus_VAJHU*p0plus_m4l_ScaleUp + bkg_prodIndep_VAMCFM*bkg_m4l_ScaleUp);
				D_bkg_prodIndep_ScaleDown = (p0plus_VAJHU*p0plus_m4l_ScaleDown) / (p0plus_VAJHU*p0plus_m4l_ScaleDown + bkg_prodIndep_VAMCFM*bkg_m4l_ScaleDown);
				D_bkg_prodIndep_ResUp     = (p0plus_VAJHU*p0plus_m4l_ResUp) / (p0plus_VAJHU*p0plus_m4l_ResUp + bkg_prodIndep_VAMCFM*bkg_m4l_ResUp);
				D_bkg_prodIndep_ResDown     = (p0plus_VAJHU*p0plus_m4l_ResDown) / (p0plus_VAJHU*p0plus_m4l_ResDown + bkg_prodIndep_VAMCFM*bkg_m4l_ResDown);

				p1plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p1plusProdIndepVA);
				p1minusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p1minusProdIndepVA);
				p1plusKD = p0plus_VAJHU / (p0plus_VAJHU + p1plus_VAJHU);
				p1minusKD = p0plus_VAJHU / (p0plus_VAJHU + p1_VAJHU);

				graviKD = p0plus_VAJHU / (p0plus_VAJHU + p2minimalVA);
				qqgraviKD = p0plus_VAJHU / (p0plus_VAJHU + p2minimalVA_qqb);
				p2h2plusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h2plusVA);
				p2h2plus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h2plusVA_qqb);
				p2h3plusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h3plusVA);
				p2h3plus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h3plusVA_qqb);
				p2hplusKD = p0plus_VAJHU / (p0plus_VAJHU + p2hplusVA);
				p2hplus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2hplusVA_qqb);
				p2bplusKD = p0plus_VAJHU / (p0plus_VAJHU + p2bplusVA);
				p2bplus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2bplusVA_qqb);
				p2h6plusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h6plusVA);
				p2h6plus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h6plusVA_qqb);
				p2h7plusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h7plusVA);
				p2h7plus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h7plusVA_qqb);
				p2hminusKD = p0plus_VAJHU / (p0plus_VAJHU + p2hminusVA);
				p2hminus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2hminusVA_qqb);
				p2h9minusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h9minusVA);
				p2h9minus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h9minusVA_qqb);
				p2h10minusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h10minusVA);
				p2h10minus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h10minusVA_qqb);
				p2mProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2mProdIndepVA);
				p2h2plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h2plusProdIndepVA);
				p2h3plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h3plusProdIndepVA);
				p2hplusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2hplusProdIndepVA);
				p2bplusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2bplusProdIndepVA);
				p2h6plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h6plusProdIndepVA);
				p2h7plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h7plusProdIndepVA);
				p2hminusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2hminusProdIndepVA);
				p2h9minusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h9minusProdIndepVA);
				p2h10minusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h10minusProdIndepVA);

			double cadd_g1g2 = 1.0;
			double cadd_g1g4 = 1.0;
			double cadd_g4ZGGG = 1.0;
			double cadd_ZG_L1 = 2.6;
			if(abs(lepIdOrdered[0])==abs(lepIdOrdered[1]) &&
				abs(lepIdOrdered[0])==abs(lepIdOrdered[2]) &&
				abs(lepIdOrdered[0])==abs(lepIdOrdered[3])){
				cadd_g1g2 = pow(1.638,2)/2.3;
				cadd_g1g4 = pow(2.521,2)/7.0;
				cadd_g4ZGGG = 1.0/7.0;
			}
			else{
				cadd_g1g2 = pow(1.638, 2) / 2.1;
				cadd_g1g4 = pow(2.521, 2) / 6.0;
				cadd_g4ZGGG = 1.0/6.0;
			};
			D_g1q2 = p0plus_VAJHU / (p0plus_VAJHU + p0_g1prime2_VAJHU);
			D_g1q2int = pg1g1prime2_VAJHU/(p0plus_VAJHU + p0_g1prime2_VAJHU);
			D_g2 = p0plus_VAJHU/(p0plus_VAJHU + p0hplus_VAJHU);
			D_g4 = p0plus_VAJHU/(p0plus_VAJHU + p0minus_VAJHU);
			D_g2int = pg1g2_VAJHU / (p0plus_VAJHU + p0hplus_VAJHU*cadd_g1g2);
			D_g4int = pg1g4_VAJHU / (p0plus_VAJHU + p0minus_VAJHU*cadd_g1g4);
			D_g2int_perp = pg1g2_pi2_VAJHU / (p0plus_VAJHU + p0hplus_VAJHU*cadd_g1g2);
			D_g4int_perp = pg1g4_pi2_VAJHU / (p0plus_VAJHU + p0minus_VAJHU*cadd_g1g4);

			D_ZG = p0plus_VAJHU / (p0plus_VAJHU + p0Zgs_VAJHU);
			D_GG = p0plus_VAJHU / (p0plus_VAJHU + p0gsgs_VAJHU);
			D_ZGint = pzzzg_VAJHU / (p0plus_VAJHU + p0Zgs_VAJHU);
			D_GGint = pzzgg_VAJHU / (p0plus_VAJHU + p0gsgs_VAJHU);
			D_ZG_PS = p0plus_VAJHU / (p0plus_VAJHU + p0Zgs_PS_VAJHU);
			D_GG_PS = p0plus_VAJHU / (p0plus_VAJHU + p0gsgs_PS_VAJHU);
			D_ZG_PSint = pzzzg_PS_VAJHU / (p0plus_VAJHU + p0Zgs_PS_VAJHU);
			D_GG_PSint = pzzgg_PS_VAJHU / (p0plus_VAJHU + p0gsgs_PS_VAJHU);

			D_ZG_L1 = p0plus_VAJHU / (p0plus_VAJHU + p0Zgs_g1prime2_VAJHU);
			D_ZG_L1int = pzzzgs_g1prime2_VAJHU*sqrt(cadd_ZG_L1) / (p0plus_VAJHU + p0Zgs_g1prime2_VAJHU*cadd_ZG_L1);
			D_ZG_L1int_perp = pzzzgs_g1prime2_pi2_VAJHU*sqrt(cadd_ZG_L1) / (p0plus_VAJHU + p0Zgs_g1prime2_VAJHU*cadd_ZG_L1);

			D_Gamma_r10=Dgg10_VAMCFM;

#ifdef LINKMELA
    if(recompute_){
      float KD, psig, pbkg;
      myMELA->computeKD(mzz,m1,m2,
		       costhetastar,
		       costheta1,
		       costheta2,
		       phi,
		       phi1,
		       KD,psig,pbkg,
		       withPt_,pt4l,
		       withY_, Y4l,
		       sqrts);
      
      melaLD = KD;
    }
#endif
    ++neventOut;
    newTree->Fill();
    if (dumpEvents) {
      if (!useJET && !VBFtag) {
	cout.setf(ios::fixed);
	int ifs=-1;
	if (outfile.Contains("4e")) ifs = 0;
	else if (outfile.Contains("4mu")) ifs = 1;
	else if (outfile.Contains("2e2mu")) ifs = 2;
	cout << run << ":" << ls << ":" << event << " " << ifs << " " << setprecision(2) << CMS_zz4l_mass << " " << m1 << " " << m2 << " " << CMS_zz4l_massErr << " "  <<  setprecision(3) << melaLD << " "  <<  setprecision(2)  << ((pt4l>0)?pt4l:-1) << " " << setprecision(3) << ((NJets_30 >= 2)?fisher:-1) << endl;
      }
    }
  }
  newTree->Write("data_obs"); 
  newFile->Close();
  
  cout << "written: " << outfile << " entries: " << neventOut << endl << endl;
  return  neventOut;
  
}


