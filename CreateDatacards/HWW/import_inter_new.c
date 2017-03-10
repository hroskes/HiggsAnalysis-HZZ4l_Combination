//#include "HZZ4L_RooSpinZeroPdf.cc+"
//#include "RooWorkSpace.h"
//#include "TFile.h"
//#include "RooRealVar.h"
using namespace RooFit;
void import_inter_new(bool jet0 = 1, bool EightTeV =1){
gROOT->SetStyle("Plain");
  //gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
  gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
  gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
  gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
  TString namebox[17]={
    "",
    "CMS_hww_MVALepEffBoundingUp",
    "CMS_hww_MVALepEffBoundingDown",
    "CMS_hww_MVALepResBoundingUp",
    "CMS_hww_MVALepResBoundingDown",
    "CMS_hww_MVAMETResBoundingUp",
    "CMS_hww_MVAMETResBoundingDown",
    "CMS_hww_MVAJESBoundingUp",
    "CMS_hww_MVAJESBoundingDown",
    "CMS_hww_VHShapeUp",
    "CMS_hww_VHShapeDown",
    "CMS_hwwof_0j_MVAggH_ALTStatBounding_8TeVUp",
    "CMS_hwwof_0j_MVAggH_ALTStatBounding_8TeVDown",
    "CMS_hwwof_0j_MVAggH_IntStatBounding_8TeVUp",
    "CMS_hwwof_0j_MVAggH_IntStatBounding_8TeVDown",
    "CMS_hwwof_0j_MVAggHStatBounding_8TeVUp",
    "CMS_hwwof_0j_MVAggHStatBounding_8TeVDown"
  };
  double T1_integral, T2_integral, T3_integral;
	double T1_statUp, T2_statUp, T3_statUp; 
	double T1_statDown, T2_statDown, T3_statDown; 
  RooWorkspace w("w");

	TString fhname = "../WWinputs/final_tree_0PM_0jets.root";
	TString fh_altname = "../WWinputs/final_tree_L1_0jets.root";
	TString fh_intername ="../WWinputs/final_tree_MixL1_n_0jets.root";
	TString fvbf_sysname="../WWinputs/SM_syst_0jet.root";

	TString fh_o= "../WWinputs/datacards/xww125p6_x125ww4l-0l1-v19/hwwof_0j.input_8TeV.root";
	TString fh_o_inter= "../WWinputs/datacards/xww125p6_x125ww4l-0l1f05ph180-v19/hwwof_0j.input_8TeV.root";

	if(!jet0){
		fhname.ReplaceAll ("0j","1j");
		fh_altname.ReplaceAll ("0j","1j");
		fh_intername.ReplaceAll ("0j","1j");
	  fh_o.ReplaceAll ("0j","1j");
		fh_o_inter.ReplaceAll ("0j","1j");
		fvbf_sysname.ReplaceAll("0j","1j");
	}
	
  TFile *fh = new TFile(fhname,"read");
  TFile *fh_inter = new TFile(fh_intername,"read");
  TFile *fh_alt = new TFile(fh_altname,"read");

	 TFile *f_o = new TFile(fh_o,"read");
	 TFile *f_o_inter = new TFile(fh_o_inter,"read");
 	TFile *fvbf_sys = new TFile(fvbf_sysname,"read");
	 
   TH1F *T1_norminal = (TH1F*)f_o->Get("histo_ggH");
   TH1F *T2_norminal = (TH1F*)f_o->Get("histo_ggH_ALT");
   TH1F *T3_norminal = (TH1F*)f_o_inter->Get("histo_ggH_ALT");

	 TH1F *ggh = (TH1F*)fvbf_sys->Get("histo1");
	 TH1F *ggh_vbf = (TH1F*)fvbf_sys->Get("histo");
  RooRealVar x ("CMS_zz4l_fg4","CMS_zz4l_fg4",0.,-1.,1.);
  RooRealVar D1 ("CMS_th1x","CMS_th1x",0, 126);

  RooArgList pdflist("pdflist"); 
  RooArgList coeffs("coeffs"); 
  TLegend *leg = new TLegend(0.65,0.53,0.85,0.83);

  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
    TH1F *T1_o = (TH1F*)fh->Get("histo");
    TH1F *T2_o = (TH1F*)fh_alt->Get("histo");
    TH1F *T3_o = (TH1F*)fh_inter->Get("histo");
	T1_integral = T1_o->Integral();
	T2_integral = T2_o->Integral();
	T3_integral = T3_o->Integral();
		double w1=0;
		double w2=0;
		double w3=0;
	for (int i =0;i<126;i++){

			w1+= T1_o->GetBinError(i+1)*T1_o->GetBinError(i+1);
			w2+= T2_o->GetBinError(i+1)*T2_o->GetBinError(i+1);
			w3+= T3_o->GetBinError(i+1)*T3_o->GetBinError(i+1);
	}
	w1= sqrt(w1);
	w2= sqrt(w2);
	w3= sqrt(w3);

  ProcessNormalization *T1_const = new ProcessNormalization(Form("T1_const_%d_%d",jet0,EightTeV),Form("T1_const_%d_%d",jet0,EightTeV),T1_integral);
  ProcessNormalization *T2_const = new ProcessNormalization(Form("T2_const_%d_%d",jet0,EightTeV),Form("T2_const_%d_%d",jet0,EightTeV),T2_integral);
  ProcessNormalization *T3_const = new ProcessNormalization(Form("T3_const_%d_%d",jet0,EightTeV),Form("T3_const_%d_%d",jet0,EightTeV),T3_integral);

	RooPlot *rplot = D1.frame();
  for (int k =0;k<11;k++){
		if (!EightTeV)
			namebox[k].ReplaceAll("8TeV","7TeV");
		if (!jet0)
			namebox[k].ReplaceAll("0j","1j");
    TString hisname = "_"+namebox[k];
    TString hisname1 = "_"+namebox[k];
    TString hisname2 = "_"+namebox[k];
    TString hisname3 = "_"+namebox[k];
    if(hisname.Contains("StatBounding")){
      if(hisname.Contains("ALT")){
				hisname1=namebox[0];
				hisname3=namebox[0];
			}
      else if (hisname.Contains("Int")){
				hisname2=namebox[0];
				hisname1=namebox[0];
				hisname3.ReplaceAll("Int","ALT");
			}
			else{
				hisname2=namebox[0];
				hisname3=namebox[0];
			}
    }
    if(k==0)
      {
				hisname1=namebox[k];
				hisname2=namebox[k];
				hisname3=namebox[k];
      }

    cout<<T1_o->Integral()<<" "<<T2_o->Integral()<< " "<<T3_o->Integral()<<endl;
    cout<<T1_norminal->Integral()<<" "<<T2_norminal->Integral()<< " "<<T3_norminal->Integral()<<endl;

    TH1F *T1= new TH1F("histo_ggH_rebin","",126,0,126);
    TH1F *T2= new TH1F("histo_ggH_ALT_rebin","",126,0,126);
    TH1F *T3= new TH1F("histo_ggH_int_rebin","",126,0,126);
    TH1F *T4= new TH1F("histo_ggH_int_test","",126,0,126);
    TH1F *T5= new TH1F("histo_ggH_int_test1","",126,0,126);
		T1->Sumw2();
		T2->Sumw2();
		T3->Sumw2();
		T4->Sumw2();
		T5->Sumw2();

		double T1_integral_o,T2_integral_o,T3_integral_o;
		
   TH1F *T1_sy = (TH1F*)f_o->Get("histo_ggH"+hisname1);
   TH1F *T2_sy = (TH1F*)f_o->Get("histo_ggH_ALT"+hisname2);
   TH1F *T3_sy = (TH1F*)f_o_inter->Get("histo_ggH_ALT"+hisname3);

    for (int i=0;i<126;i++){
			double r1=1;
			double r2=1;
			double r3=1;
			if (!hisname.Contains("StatBounding")){
		if (hisname.Contains("VH")){
			if (ggh->GetBinContent(i+1)!=0){
			if (hisname.Contains("Up")){
			r1 = ggh_vbf->GetBinContent(i+1)/ggh->GetBinContent(i+1); 
			r2 = ggh_vbf->GetBinContent(i+1)/ggh->GetBinContent(i+1); 
			r3 = ggh_vbf->GetBinContent(i+1)/ggh->GetBinContent(i+1); 
			}
			else{
			r1 = 2-ggh_vbf->GetBinContent(i+1)/ggh->GetBinContent(i+1); 
			r2 = 2-ggh_vbf->GetBinContent(i+1)/ggh->GetBinContent(i+1); 
			r3 = 2-ggh_vbf->GetBinContent(i+1)/ggh->GetBinContent(i+1); 
			if (r1<0){
				r1=0; r2=0; r3=0;
			}
			}
			}
		}
		
		else if(T1_norminal->GetBinContent(i+1)!=0 && T2_norminal->GetBinContent(i+1)!=0 && T3_norminal->GetBinContent(i+1)!=0){
//		if(T1_norminal->GetBinContent(i+1)!=0 && T2_norminal->GetBinContent(i+1)!=0 && T3_norminal->GetBinContent(i+1)!=0){
			r1= T1_sy->GetBinContent(i+1)/T1_norminal->GetBinContent(i+1);
			r2= T2_sy->GetBinContent(i+1)/T2_norminal->GetBinContent(i+1);
			r3= T3_sy->GetBinContent(i+1)/T3_norminal->GetBinContent(i+1);
		}

      T1->SetBinContent(i+1,T1_o->GetBinContent(i+1)*r1);
      T2->SetBinContent(i+1,T2_o->GetBinContent(i+1)*r2);
 	    T3->SetBinContent(i+1,T3_o->GetBinContent(i+1)*r3);

      T1->SetBinError(i+1,T1_o->GetBinError(i+1)*r1);
      T2->SetBinError(i+1,T2_o->GetBinError(i+1)*r2);
 	    T3->SetBinError(i+1,T3_o->GetBinError(i+1)*r3);
			}

			else {
				T1->SetBinContent(i+1,T1_o->GetBinContent(i+1));
				T2->SetBinContent(i+1,T2_o->GetBinContent(i+1));
				T3->SetBinContent(i+1,T3_o->GetBinContent(i+1));
			  if (hisname.Contains("Int")){
					if (hisname.Contains("Up"))
						T3->SetBinContent(i+1,T3_o->GetBinContent(i+1) * (1+w3/T3_integral ));
					else 
						T3->SetBinContent(i+1,T3_o->GetBinContent(i+1) * (1-w3/T3_integral ) );
				}
			  else if (hisname.Contains("ALT")){
					if (hisname.Contains("Up"))
						T2->SetBinContent(i+1,T2_o->GetBinContent(i+1) * (1+w2/T2_integral) );
					else 
						T2->SetBinContent(i+1,T2_o->GetBinContent(i+1) * (1-w2/T2_integral) );
				}
			  else {
					if (hisname.Contains("Up"))
						T1->SetBinContent(i+1,T1_o->GetBinContent(i+1) * (1+w1/T1_integral) );
					else 
						T1->SetBinContent(i+1,T1_o->GetBinContent(i+1) * (1-w1/T1_integral) );
				}
		//		cout<< T1_o->GetBinError(i+1)<< " "<< T1_norminal->GetBinError(i+1)<<endl;
				//cout<< T1->GetBinContent(i+1)<<" "<< T1_o->GetBinContent(i+1)<<endl;
				//cout<< T2->GetBinContent(i+1)<<" "<< T2_o->GetBinContent(i+1)<<endl;
				//cout<< T3->GetBinContent(i+1)<<" "<< T3_o->GetBinContent(i+1)<<endl;
			}

//			cout<<r1<<" "<<r2<<" "<<r3<<endl;
      T1->SetBinError(i+1,T1_o->GetBinError(i+1)*r1);
      T2->SetBinError(i+1,T2_o->GetBinError(i+1)*r2);
 	    T3->SetBinError(i+1,T3_o->GetBinError(i+1)*r3);
			}
			 
    cout<<T1->Integral()<<" "<<T2->Integral()<< " "<<T3->Integral()<<endl;
		if (hisname.Contains("Up")){
			T1_statUp = T1->Integral(); 
			T2_statUp = T2->Integral(); 
			T3_statUp = T3->Integral(); 
		}
		else if (hisname.Contains("Down")){
	      T1_statDown = T1->Integral();
	      T2_statDown = T2->Integral();
	      T3_statDown = T3->Integral();
		}
			 T3->Add(T1,-1);
			 T3->Add(T2,-1);
//			 T3->Scale(-1);
			 
		if(k==0){
			T1->Draw("histe");
			T1->SetMaximum(T1->GetMaximum()*1.3);
			T1->SetLineColor(1);
			T2->SetLineColor(2);
			T1->SetMarkerColor(1);
			T2->SetMarkerColor(2);
			T2->Draw("histesame");
			leg->AddEntry(T1,"0m+,ggH");
			leg->AddEntry(T2,"alt,ggH");
		}
			if (k==9){
				T1->SetLineColor(8);
				T2->SetLineColor(4);
				T1->SetMarkerColor(8);
				T2->SetMarkerColor(4);
				T1->Draw("histesame");
				T2->Draw("histesame");
			leg->AddEntry(T1,"0m+,ggH+all");
			leg->AddEntry(T2,"alt,ggH+all");
			leg->Draw();
			if(jet0){
				gPad->Print("flambda1_template_0jet.eps");
				gPad->Print("flambda1_template_0jet.png");
			}
			else{
				gPad->Print("flambda1_template_1jet.eps");
				gPad->Print("flambda1_template_1jet.png");
			}
			}
			//T3->Draw("same");
			//T1->SetMinimum(T3->GetMinimum());
		
			//T4->Reset();
			//T4->Add(T1,0.5);
			//T4->Add(T2,0.5);
			//T4->Add(T3,-0.5);
			//T4->SetLineColor(5);
			//T4->Draw("same");

			//T5->Reset();
			//T5->Add(T1,0.6);
			//T5->Add(T2,0.4);
			//T5->Add(T3,0.4899);
			//T5->SetLineColor(8);
			//T5->Draw("same");
			//gPad->Print("template.png");
			//for (int i =0;i<126;i++){
			//if (T5->GetBinContent(i+1)<0)
			//		cout<<i<<" "<<T5->GetBinContent(i+1)<<" "<<T5->GetBinError(i+1)<<endl;
			//}

		//}
		
				
    float dLowX = T1->GetXaxis()->GetXmin();
    float dHighX = T1->GetXaxis()->GetXmax();
    RooDataHist *hist1 = new RooDataHist("hist1"+hisname,"hist1"+hisname,RooArgList(D1),T1);
    RooDataHist *hist2 = new RooDataHist("hist2"+hisname,"hist2"+hisname,RooArgList(D1),T2);
    RooDataHist *hist3 = new RooDataHist("hist3"+hisname,"hist3"+hisname,RooArgList(D1),T3);
    RooHistFunc *histfunc1 = new RooHistFunc("histfunc1"+hisname,"histfunct1"+hisname,RooArgSet(D1),*hist1);
    RooHistFunc *histfunc2 = new RooHistFunc("histfunc2"+hisname,"histfunct2"+hisname,RooArgSet(D1),*hist2);
    RooHistFunc *histfunc3 = new RooHistFunc("histfunc3"+hisname,"histfunct3"+hisname,RooArgSet(D1),*hist3);
    HZZ4L_RooSpinZeroPdf_1D *ggHpdf = new HZZ4L_RooSpinZeroPdf_1D("tmpggH"+hisname,"tmpggH"+hisname,D1,x,RooArgList(*histfunc1,*histfunc2,*histfunc3));
//    w.import(*ggHpdf,RecycleConflictNodes());
		pdflist.add(*ggHpdf);

	if(k!=0 && k%2==0 ){
		TString sysname = namebox[k];
		sysname. ReplaceAll("Down","");
		cout<<sysname<<endl;
		RooRealVar *sysrr=new RooRealVar(sysname,sysname,-7,7);	
		coeffs.add(*sysrr);
//		cout<<"Gaussian::"+sysname+"_Pdf("+sysname+"[-7,7], "+sysname+"_In[0,-7,7], 1)"<<endl;
//		w.factory("Gaussian::"+sysname+"_Pdf("+sysname+"[-7,7], "+sysname+"_In[0,-7,7], 1)");

		cout<<T1_statDown/T1_integral<<" "<<T1_statUp/T1_integral<<endl;
		cout<<T2_statDown/T2_integral<<" "<<T2_statUp/T2_integral<<endl;
		cout<<T3_statDown/T3_integral<<" "<<T3_statUp/T3_integral<<endl;
		if(k!=10){
		T1_const->addAsymmLogNormal(T1_statDown/T1_integral, T1_statUp/T1_integral, *sysrr);	
		T2_const->addAsymmLogNormal(T2_statDown/T2_integral, T2_statUp/T2_integral, *sysrr);	
		T3_const->addAsymmLogNormal(T3_statDown/T3_integral, T3_statUp/T3_integral, *sysrr);	
		}
	}

  }


  RooFormulaVar *ggH_norm=new RooFormulaVar("ggH__norm","((1-abs(@0))*@1 + abs(@0)*@2 + (@0>0 ? 1.:-1.)*sqrt(abs(@0)*(1-abs(@0)))*(@3-@1-@2))/@1 >0 ? ((1-abs(@0))*@1 + abs(@0)*@2 + (@0>0 ? 1.:-1.)*sqrt(abs(@0)*(1-abs(@0)))*(@3-@1-@2))/@1 : 1.e-20", RooArgList(x, *T1_const, *T2_const, *T3_const));
pdflist.Print("v");	
coeffs.Print("v");	
	VerticalInterpPdf *ggH = new VerticalInterpPdf("ggH_","",pdflist,coeffs,1,0);
    w.import(*ggH,RecycleConflictNodes());
	 w.import(*ggH_norm,RecycleConflictNodes());
  w.var("CMS_zz4l_fg4").setVal(0);
	
	TString jetName= jet0 ? "0j_":"1j_";
	TString EnergyName = EightTeV ? "_8TeV":"7TeV";
	TString outfileName = "signal_inter_"+ jetName + EnergyName + ".root"; 
  w.writeToFile(outfileName);
}
