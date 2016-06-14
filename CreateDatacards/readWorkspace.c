
//#include "HZZ2L2QRooPdfs.cc+"
//#include "HZZ4L_RooHighmass.cc+"
//#include "RooHighmass_conv.cc+"
//#include "ProcessNormalization.cc+"
using namespace RooFit;
using namespace RooStats;
void readWorkspace(TString file){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
  gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
  gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
//	TFile *f=new TFile("/scratch0/hep/hroskes/Summer2015_VBF/forCJLST/CMSSW_6_1_1/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards/cards_newyields_newsignal_oldbkg/HCG/125/fixedMu_8TeV.root");
	TFile *f=new TFile(file);
	RooWorkspace *w= (RooWorkspace*) f->Get("w");
	w->Print();
  w->exportToCint("w");
RooPlot *frame = w::CMS_zz4l_fg4->frame();
RooPlot *frame_kd = w::D_0minus_decay->frame();
//w::shapeSig_ggH_ch1__norm->plotOn(frame);
//frame->Draw();
w::CMS_zz4l_fg4->setVal(0);
w::ggH->plotOn(frame_kd);
w::CMS_zz4l_fg4->setVal(0.5);
w::ggH->plotOn(frame_kd,LineColor(2));
w::CMS_zz4l_fg4->setVal(-0.5);
w::ggH->plotOn(frame_kd,LineColor(8));
w::CMS_zz4l_fg4->setVal(1);
w::ggH->plotOn(frame_kd,LineColor(1));
frame_kd->Draw();
gPad->Print("orig_tempaltes.png");

}
