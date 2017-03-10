#include <iostream>
#include <cmath>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


void modifySigTemplates(TString dir, TString sqrts = "7TeV", TString channame = "2e2mu"){
	TString cinputNominal = dir + "/" + sqrts + "/" + channame + "_templates_Nominal.root";
	TString cinputSysUp = dir + "/" + sqrts + "/" + channame + "_templates_SysUp.root";
	TString cinputSysDown = dir + "/" + sqrts + "/" + channame + "_templates_SysDown.root";
	TString cinput[3] = {
		cinputNominal,
		cinputSysUp,
		cinputSysDown
	};

	TString coutput = dir + "/" + sqrts + "/" + channame + "_templates_Modified_Nominal_ScaleResUpDown.root";
	TFile* foutput = new TFile(coutput, "recreate");

	const int kNumTemplates=9;
	const int kNumSysts=3;

	TFile* finput[kNumSysts];
	TH2F* hsig[kNumSysts * kNumTemplates];
	TH3F* hsig_scaleres[kNumSysts * kNumTemplates];
	for (int f = 0; f < kNumSysts; f++){
		finput[f] = new TFile(cinput[f], "read");
		cout << "Reading file: " << finput[f]->GetName() << "..." << endl;

		string ctemplate_main = "T_2D_";
		string ctemplate_main_new = "T_3D_";
		for (int t = 0; t < kNumTemplates; t++){
			char tcode[2];
			sprintf(tcode, "%i", t+1);
			string ctemplate = ctemplate_main + tcode;
			string ctemplate_new = ctemplate_main_new + tcode;
			string ctemplate_scaleres = ctemplate_new + "_ScaleRes";
			if ((f - 1) % 2 == 0) ctemplate_scaleres = ctemplate_scaleres + "Up";
			if ((f - 1) % 2 == 1) ctemplate_scaleres = ctemplate_scaleres + "Down";
			hsig[kNumTemplates * f + t] = (TH2F*) finput[f]->Get(ctemplate.c_str());
//			if(t==2) hsig[kNumTemplates * f + t]->Scale(1./10.);

			if (f>0){
				if ( (hsig[kNumTemplates * f + t]->Integral()) != 0) hsig[kNumTemplates * f + t]->Scale( fabs(hsig[kNumTemplates * 0 + t]->Integral()) / fabs(hsig[kNumTemplates * f + t]->Integral()) );
			}
			int nbinsX = hsig[kNumTemplates * f + t]->GetNbinsX();
			int nbinsY = hsig[kNumTemplates * f + t]->GetNbinsY();
			double xmin = hsig[kNumTemplates * f + t]->GetXaxis()->GetXmin();
			double xmax = hsig[kNumTemplates * f + t]->GetXaxis()->GetXmax();
			double ymin = hsig[kNumTemplates * f + t]->GetYaxis()->GetXmin();
			double ymax = hsig[kNumTemplates * f + t]->GetYaxis()->GetXmax();

			if(f>0) hsig_scaleres[kNumTemplates * f + t] = new TH3F(
				ctemplate_scaleres.c_str(),
				ctemplate_scaleres.c_str(),
				nbinsX,xmin,xmax,
				nbinsY,ymin,ymax,
				1,0,1
				);
			else hsig_scaleres[kNumTemplates * f + t] = new TH3F(
				ctemplate_new.c_str(),
				ctemplate_new.c_str(),
				nbinsX,xmin,xmax,
				nbinsY,ymin,ymax,
				1,0,1
				);
			for (int binx = 1; binx <= nbinsX; binx++){
				for (int biny = 1; biny <= nbinsY; biny++){
					double bincontent = hsig[kNumTemplates * f + t]->GetBinContent(binx, biny);
					if(bincontent!=bincontent) cout << "Bin content is NaN! Template: " << tcode << endl;
					for (int binz = 1; binz <= 1; binz++) hsig_scaleres[kNumTemplates * f + t]->SetBinContent(binx, biny, binz, bincontent);
				}
			}
		};
	};
	cout << "Integrals:\n"
		<< "Template\tNominal\tUp\tDown" << endl;
	for (int t = 0; t < kNumTemplates; t++){
		TString ccanvas = "validateSig_";
		char tcode[2];
		sprintf(tcode, "%i", t+1);
		ccanvas = ccanvas + tcode;

		cout << t+1 << '\t';
		for (int f = 0; f < kNumSysts; f++){
			foutput->WriteTObject(hsig_scaleres[kNumTemplates * f + t]);
			cout << hsig_scaleres[kNumTemplates * f + t]->Integral() << '\t';
		};
		cout << endl;

		TH2F* h2dproj[kNumSysts];

		TString ccanvas_xy = ccanvas + "_xy";
		TCanvas* c = new TCanvas(ccanvas_xy,"",1000,800);
		c->Divide(2,2);
		double Ncontrib = hsig_scaleres[kNumTemplates * 0 + t]->Integral();
		for (int f = 0; f < kNumSysts; f++){
			c->cd(f + 1);
			h2dproj[f] = new TH2F(Form("hproj_xy_%i",f),"",
				hsig_scaleres[kNumTemplates * f + t]->GetNbinsX(),hsig_scaleres[kNumTemplates * f + t]->GetXaxis()->GetXmin(),hsig_scaleres[kNumTemplates * f + t]->GetXaxis()->GetXmax(),
				hsig_scaleres[kNumTemplates * f + t]->GetNbinsY(),hsig_scaleres[kNumTemplates * f + t]->GetYaxis()->GetXmin(),hsig_scaleres[kNumTemplates * f + t]->GetYaxis()->GetXmax()
				);
			for (int binx = 1; binx <= hsig_scaleres[kNumTemplates * f + t]->GetNbinsX(); binx++){
				for (int biny = 1; biny <= hsig_scaleres[kNumTemplates * f + t]->GetNbinsY(); biny++){
					double sum = 0;
					for (int binz = 1; binz <= hsig_scaleres[kNumTemplates * f + t]->GetNbinsZ(); binz++) sum += hsig_scaleres[kNumTemplates * f + t]->GetBinContent(binx, biny, binz);
					h2dproj[f]->SetBinContent(binx,biny,sum);
				}
			}
			h2dproj[f]->GetXaxis()->SetTitle("m_{4l} (GeV)");
			h2dproj[f]->GetYaxis()->SetTitle("#it{D}^{kin}_{bkg}");
			h2dproj[f]->Scale( 1.0/fabs( h2dproj[f]->Integral() ) );
			if(f==0) h2dproj[f]->SetTitle(Form("T %s Nominal",tcode));
			if(f==1) h2dproj[f]->SetTitle(Form("T %s Up",tcode));
			if(f==2) h2dproj[f]->SetTitle(Form("T %s Down",tcode));
			h2dproj[f]->Draw("colz");
		}
		foutput->WriteTObject(c);
		delete c;
		for (int f = 0; f < kNumSysts; f++){
			delete h2dproj[f];
		}

		TString ccanvas_xz = ccanvas + "_zx";
		c = new TCanvas(ccanvas_xz,"",1000,800);
		c->Divide(2,2);
		c->cd(1);
		hsig_scaleres[kNumTemplates * 0 + t]->Project3D("xz")->Draw("colz");
		c->cd(2);
		hsig_scaleres[kNumTemplates * 1 + t]->Project3D("xz")->Draw("colz");
		c->cd(3);
		hsig_scaleres[kNumTemplates * 2 + t]->Project3D("xz")->Draw("colz");
		foutput->WriteTObject(c);
		delete c;
		TString ccanvas_yz = ccanvas + "_yz";
		c = new TCanvas(ccanvas_yz,"",1000,800);
		c->Divide(2,2);
		c->cd(1);
		hsig_scaleres[kNumTemplates * 0 + t]->Project3D("zy")->Draw("colz");
		c->cd(2);
		hsig_scaleres[kNumTemplates * 1 + t]->Project3D("zy")->Draw("colz");
		c->cd(3);
		hsig_scaleres[kNumTemplates * 2 + t]->Project3D("zy")->Draw("colz");
		foutput->WriteTObject(c);
		delete c;
	};

	foutput->Close();
	for (int f = 0; f < kNumSysts; f++) finput[f]->Close();
}


void modifyBkgTemplates (TString dir, TString sqrts = "7TeV", TString channame = "2e2mu"){
	TString cinput = dir + "/" + sqrts + "/" + channame + "_templates_Nominal_bkg.root";
	TString coutput = dir + "/" + sqrts + "/" + channame + "_templates_Modified_Nominal_bkg.root";
	
	TFile* finput = new TFile(cinput,"read");
	TFile* foutput = new TFile(coutput,"recreate");
	TH2F* ggzz_2d = (TH2F*) finput->Get("template_ggZZ");
	ggzz_2d->SetName("template_ggZZ_2D");
	TH2F* zx_2d= (TH2F*) finput->Get("template_ZX");
	zx_2d->SetName("template_ZX_2D");

	TH3F* ggzz = new TH3F("template_ggZZ","",
		ggzz_2d->GetNbinsX(),ggzz_2d->GetXaxis()->GetXmin(),ggzz_2d->GetXaxis()->GetXmax(),
		ggzz_2d->GetNbinsY(),ggzz_2d->GetYaxis()->GetXmin(),ggzz_2d->GetYaxis()->GetXmax(),
		1,0,1);
	TH3F* zx = new TH3F("template_ZX","",
		zx_2d->GetNbinsX(),zx_2d->GetXaxis()->GetXmin(),zx_2d->GetXaxis()->GetXmax(),
		zx_2d->GetNbinsY(),zx_2d->GetYaxis()->GetXmin(),zx_2d->GetYaxis()->GetXmax(),
		1,0,1);
	for (int binz = 1; binz <= ggzz->GetNbinsZ(); binz++){
		for (int binx = 1; binx <= ggzz->GetNbinsX(); binx++){
			for (int biny = 1; biny <= ggzz->GetNbinsY(); biny++){
				double bincontent = ggzz_2d->GetBinContent(binx, biny);
				ggzz->SetBinContent(binx, biny, binz, bincontent);
			}
		}
	}
	for (int binz = 1; binz <= zx->GetNbinsZ(); binz++){
		for (int binx = 1; binx <= zx->GetNbinsX(); binx++){
			for (int biny = 1; biny <= zx->GetNbinsY(); biny++){
				double bincontent = zx_2d->GetBinContent(binx, biny);
				zx->SetBinContent(binx, biny, binz, bincontent);
			}
		}
	}
	ggzz->Scale(1.0/ggzz->Integral());
	zx->Scale(1.0/zx->Integral());

	TH3F* zx_up = (TH3F*) ggzz->Clone("template_ZX_Up");
	TH3F* zx_down = (TH3F*) zx->Clone("template_ZX_Down");
	zx_up->SetNameTitle("template_ZX_Up","template_ZX_Up");
	zx_down->SetNameTitle("template_ZX_Down","template_ZX_Down");
	for (int binx = 1; binx <= zx_down->GetNbinsX(); binx++){
		for (int biny = 1; biny <= zx_down->GetNbinsY(); biny++){
			for (int binz = 1; binz <= zx_down->GetNbinsZ(); binz++){
				double bincontent_zx = zx->GetBinContent(binx, biny, binz);
				double bincontent_zxup = zx_up->GetBinContent(binx, biny, binz);
				double bincontent = 2.0*bincontent_zx - bincontent_zxup;
				if (bincontent <= 0) bincontent = 1.0e-20;
				zx_down->SetBinContent(binx, biny, binz,bincontent);
			};
		};
	};
	zx_down->Scale(1.0/zx_down->Integral());

	foutput->WriteTObject(ggzz);
	foutput->WriteTObject(zx);
	foutput->WriteTObject(zx_up);
	foutput->WriteTObject(zx_down);

	TString ccanvas = "validateBkg";
	TString ccanvas_xy = ccanvas + "_xy";
	TCanvas* c = new TCanvas(ccanvas_xy,"",1000,800);
	c->Divide(2,2);
	c->cd(2);
	gStyle->SetPalette(1);
	zx->SetLineColor(1);
	zx->Project3D("yx")->Draw("colz");
	c->cd(1);
	ggzz->SetLineColor(2);
	ggzz->Project3D("yx")->Draw("colz");
	c->cd(3);
	zx_up->SetLineColor(4);
	zx_up->Project3D("yx")->Draw("colz");
	c->cd(4);
	zx_down->SetLineColor(4);
	zx_down->Project3D("yx")->Draw("colz");
	foutput->WriteTObject(c);
	delete c;
	TString ccanvas_xz = ccanvas + "_zx";
	c = new TCanvas(ccanvas_xz,"",1000,800);
	c->Divide(2,2);
	c->cd(2);
	gStyle->SetPalette(1);
	zx->SetLineColor(1);
	zx->Project3D("xz")->Draw("colz");
	c->cd(1);
	ggzz->SetLineColor(2);
	ggzz->Project3D("xz")->Draw("colz");
	c->cd(3);
	zx_up->SetLineColor(4);
	zx_up->Project3D("xz")->Draw("colz");
	c->cd(4);
	zx_down->SetLineColor(4);
	zx_down->Project3D("xz")->Draw("colz");
	foutput->WriteTObject(c);
	delete c;
	TString ccanvas_yz = ccanvas + "_yz";
	c = new TCanvas(ccanvas_yz,"",1000,800);
	c->Divide(2,2);
	c->cd(2);
	gStyle->SetPalette(1);
	zx->SetLineColor(1);
	zx->Project3D("zy")->Draw("colz");
	c->cd(1);
	ggzz->SetLineColor(2);
	ggzz->Project3D("zy")->Draw("colz");
	c->cd(3);
	zx_up->SetLineColor(4);
	zx_up->Project3D("zy")->Draw("colz");
	c->cd(4);
	zx_down->SetLineColor(4);
	zx_down->Project3D("zy")->Draw("colz");
	foutput->WriteTObject(c);
	delete c;

	foutput->Close();
	finput->Close();
}

void modifyTemplates_Z4lBuilder(TString dir, int processSig=1, int processBkg=1){
	gStyle->SetPadRightMargin(0.20);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadLeftMargin(0.16);
	gROOT->ForceStyle();

	TString channame[3] = { "4mu","4e","2e2mu" };
	TString sqrts[2] = { "7TeV","8TeV" };
	for (int s = 0; s < 2; s++){
		for (int f = 0; f < 3; f++){
			if(processBkg==1) modifyBkgTemplates(dir, sqrts[s], channame[f]);
			if(processSig==1) modifySigTemplates(dir, sqrts[s], channame[f]);
		};
	};
}