#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TColor.h"
#include "TString.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

using namespace std;

void countTauContamination(int erg_tev=7){
	string user_dir = "/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/";
	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV/",erg_tev);
	char* user_folder[3]={
		"4mu",
		"4e",
		"2mu2e"
	};
	char* cQQBZZ[6] = {
		"HZZ4lTree_ZZTo4mu.root",
		"HZZ4lTree_ZZTo4e.root",
		"HZZ4lTree_ZZTo2e2mu.root",
		"HZZ4lTree_ZZTo2mu2tau.root",
		"HZZ4lTree_ZZTo2e2tau.root",
		"HZZ4lTree_ZZTo4tau.root"
	};

	cout << "Channel" << '\t';
	cout << "xsec" << '\t'
		<< "ngen_all" << '\t'
		<< "ngen" << '\t'
		<< "nreco_all" << '\t'
		<< "nreco" << endl;

	for (int f = 0; f < 3; f++){
		cout << user_folder[f] << endl;
		for (int d = 0; d < 6; d++){
			char cinput[1000];
			sprintf(cinput, "%s%s%s%s%s", user_dir.c_str(),erg_dir, user_folder[f],"/", cQQBZZ[d]);
			TFile* finput = new TFile(cinput, "read");
			TTree* t = (TTree*)finput->Get("SelectedTree");
			TH1F* h1D = (TH1F*)finput->Get("hCounters");
			TH2F* h2D = (TH2F*)finput->Get("Counter_ZZQQB_STU");
			TH1F* h = new TH1F("h", "", 1, 80, 100);
			TH1F* hxsec = new TH1F("hxsec", "", 1, 80, 100);
			TH1F* hall = new TH1F("hall", "", 1, 0, 3000);
			t->Draw("ZZMass>>h", "MC_weight_noxsec");
			t->Draw("ZZMass>>hall", "MC_weight_noxsec");
			t->Draw("ZZMass>>hxsec", "MC_weight");
			double hint = h->Integral();
			double hallint = hall->Integral();
			double hxsecint = hxsec->Integral();
			double ngen = h2D->Integral(2, 3, 1, 1)/h1D->GetBinContent(1);
			double ngenall = h2D->Integral(1, h2D->GetNbinsX(), 1, 1)/h1D->GetBinContent(1);
			double xsec = (hint!=0 ? hxsecint / hint : 0);

			cout << cQQBZZ[d] << '\t';
			cout << xsec << '\t'
				<< ngenall << '\t'
				<< ngen << '\t'
				<< hallint << '\t'
				<< hint << endl;

			delete hall;
			delete hxsec;
			delete h;
			delete h2D;
			delete h1D;
			delete t;
			finput->Close();
		}
	}
}