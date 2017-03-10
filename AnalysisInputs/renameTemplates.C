#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include "TSystem.h"


using namespace std;

const TString cinput_main = "/work-zfs/lhc/usarica/CMS-related/SpinParityTemplates_RunI/";
const TString coutput_main = "../CreateDatacards/templates2D/";

void renameTemplates_sig(TString infile, TString outdir, TString outfile, bool input2D=false, bool output2D=true);
void renameTemplates_bkg(TString infile, TString outdir, TString outfile);

void renameTemplates(){
  TString channame[3]={
    "2e2mu", "4e", "4mu"
  };
  for (unsigned int CoM=7; CoM<=8; CoM++){
    for (unsigned int ichan=0; ichan<3; ichan++){
      renameTemplates_sig(Form("fa2phaseAdap_new_%iTeV/", CoM) + channame[ichan] + "_fa2phaseAdap_new", Form("fa2_3D_YR3/%iTeV/", CoM), channame[ichan] + "_templates_Modified", false, true);
      renameTemplates_sig(Form("fa3Adap_new_%iTeV/", CoM) + channame[ichan] + "_fa3Adap_new", Form("fa3_3D_YR3/%iTeV/", CoM), channame[ichan] + "_templates_Modified", false, true);

      renameTemplates_bkg(Form("fa2phaseAdap_new_%iTeV/", CoM) + channame[ichan] + "_fa2phaseAdap_new", Form("fa2_3D_YR3/%iTeV/", CoM), channame[ichan] + "_templates_Modified");
      renameTemplates_bkg(Form("fa3Adap_new_%iTeV/", CoM) + channame[ichan] + "_fa3Adap_new", Form("fa3_3D_YR3/%iTeV/", CoM), channame[ichan] + "_templates_Modified");
    }
  }
}


void renameTemplates_sig(TString infile, TString outdir, TString outfile, bool input2D, bool output2D){
  TString cinput = cinput_main + infile;
  TString coutput = coutput_main + outdir + outfile;
  gSystem->Exec("mkdir -p " + coutput_main + outdir);

  TFile* fout_sig = TFile::Open(coutput+"_Nominal_ScaleResUpDown.root", "recreate");
  TFile* fin_sig[3] ={
    TFile::Open(cinput+".root", "read"),
    TFile::Open(cinput+"_ResScaleUp.root", "read"),
    TFile::Open(cinput+"_ResScaleDown.root", "read")
  };

  vector<TString> cin_sig;
  vector<TString> cout_sig;
  if (!input2D){
    cin_sig.push_back("template0PlusAdapSmoothMirror");
    cin_sig.push_back("template0MinusAdapSmoothMirror");
    cin_sig.push_back("templateIntAdapSmoothMirror");
    cin_sig.push_back("templateIntPi2AdapSmoothMirror");
  }
  else{
    cin_sig.push_back("template0PlusAdapSmoothMirror");
    cin_sig.push_back("template0HPlusAdapSmoothMirror");
    cin_sig.push_back("templateInt1AdapSmoothMirror");
    cin_sig.push_back("templateInt1Pi2AdapSmoothMirror");

    cin_sig.push_back("template0MinusAdapSmoothMirror");
    cin_sig.push_back("templateInt2AdapSmoothMirror");
    cin_sig.push_back("templateInt2Pi2AdapSmoothMirror");
    cin_sig.push_back("templateInt3AdapSmoothMirror");
    cin_sig.push_back("templateInt3Pi2AdapSmoothMirror");
  }
  if (output2D){
    cout_sig.push_back("T_3D_1");
    cout_sig.push_back("T_3D_2");
    cout_sig.push_back("T_3D_4"); // 1+2
    cout_sig.push_back("T_3D_7"); // 1+2

    cout_sig.push_back("T_3D_3");
    cout_sig.push_back("T_3D_5"); // 1+3
    cout_sig.push_back("T_3D_8"); // 1+3
    cout_sig.push_back("T_3D_6"); // 2+3
    cout_sig.push_back("T_3D_9"); // 2+3
  }
  else{
    cout_sig.push_back("T_3D_1");
    cout_sig.push_back("T_3D_2");
    cout_sig.push_back("T_3D_4");
  }

  vector<TH3F*> hinput;
  for (unsigned int f=0; f<3; f++){
    for (unsigned int ih=0; ih<cout_sig.size(); ih++){
      TH3F* hin = 0;
      if (ih<cin_sig.size()) hin = (TH3F*)fin_sig[f]->Get(cin_sig.at(ih));
      if (hin==0){
        if (ih<cin_sig.size()) cout << "File " << fin_sig[f]->GetName() << " does not contain histogram " << cin_sig[ih] << endl;
        else cout << "File " << fin_sig[f]->GetName() << " would not contain a histogram that corresponds to " << cout_sig[ih] << endl;
      }
      if (hin==0 && ih>0){
        hin = (TH3F*)hinput.at(0)->Clone(cout_sig.at(ih));
        hin->Reset("ICES");
      }
      else if (hin==0) assert(0);
      if (f==0) hin->SetName(cout_sig.at(ih));
      else if (f==1) hin->SetName(cout_sig.at(ih) + "_ScaleResUp");
      else hin->SetName(cout_sig.at(ih) + "_ScaleResDown");
      hinput.push_back(hin);
    }
  }
  for (unsigned int ih=0; ih<hinput.size(); ih++){
    cout << "Final integral of " << hinput.at(ih)->GetName() << " = " << hinput.at(ih)->Integral() << endl;
    fout_sig->WriteTObject(hinput.at(ih));
  }
  for (unsigned int ih=0; ih<hinput.size(); ih++) delete hinput.at(ih);
  hinput.clear();
  for (unsigned int f=0; f<3; f++) fin_sig[f]->Close();
  fout_sig->Close();
}

void renameTemplates_bkg(TString infile, TString outdir, TString outfile){
  TString cinput = cinput_main + infile;
  TString coutput = coutput_main + outdir + outfile;
  gSystem->Exec("mkdir -p " + coutput_main + outdir);

  TFile* fout_bkg = TFile::Open(coutput+"_Nominal_bkg.root", "recreate");
  TFile* fin_bkg = TFile::Open(cinput+"_bkg.root", "read");

  vector<TString> cin_bkg;
  vector<TString> cout_bkg;

  cin_bkg.push_back("template_ZX");
  cin_bkg.push_back("template_qqZZ");
  cin_bkg.push_back("template_ggZZ");

  cout_bkg.push_back("template_ZX");
  cout_bkg.push_back("template_qqZZ");
  cout_bkg.push_back("template_ggZZ");
  cout_bkg.push_back("template_ZX_Down");
  cout_bkg.push_back("template_ZX_Up");

  vector<TH3F*> hinput;
  for (unsigned int ih=0; ih<cout_bkg.size(); ih++){
    TH3F* hin = 0;
    if (ih<cin_bkg.size()){
      hin = (TH3F*)fin_bkg->Get(cin_bkg.at(ih));
      if (hin==0){ cout << "File " << fin_bkg->GetName() << " does not contain histogram " << cin_bkg[ih] << endl; assert(0); }
    }
    else{
      hin = (TH3F*)hinput.at(0)->Clone(cout_bkg.at(ih)); hin->Add(hinput.at(ih-3), 1.); hin->Add(hinput.at(4-ih), -1.);
      double pre_integral = hin->Integral();
      for (int ix=1; ix<hin->GetNbinsX()+1; ix++){
        for (int iy=1; iy<hin->GetNbinsY()+1; iy++){
          for (int iz=1; iz<hin->GetNbinsZ()+1; iz++){
            double bincontent = hin->GetBinContent(ix, iy, iz);
            if (bincontent<=0.){
              //cout << "Bin (" << ix << "," << iy << "," << iz << ") of histogram " << cout_bkg.at(ih) << " is negative(" << bincontent << "). Fixing with value ";
              bincontent = pre_integral/((double)hin->GetNbinsX()*hin->GetNbinsY()*hin->GetNbinsZ())*1e-15;
              hin->SetBinContent(ix, iy, iz, bincontent);
              //cout << bincontent << endl;
            }
            hin->SetBinError(ix, iy, iz, 0);
          }
        }
      }
      double post_integral = hin->Integral();
      cout << "Scaling " << cout_bkg.at(ih) << " by " << pre_integral  << " / " << post_integral << endl;
      hin->Scale(pre_integral/post_integral);
    }
    hin->SetName(cout_bkg.at(ih));
    hinput.push_back(hin);
  }
  for (unsigned int ih=0; ih<hinput.size(); ih++){
    cout << "Final integral of " << hinput.at(ih)->GetName() << " = " << hinput.at(ih)->Integral() << endl;
    fout_bkg->WriteTObject(hinput.at(ih));
  }
  for (unsigned int ih=0; ih<hinput.size(); ih++) delete hinput.at(ih);
  hinput.clear();
  fin_bkg->Close();
  fout_bkg->Close();
}
