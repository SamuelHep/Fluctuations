#include <string>
#include <iostream>

using namespace std;

void GlauberFit(TString glauberFileName, TString multFileName, 
		TString histname, TString chisqmethod, 
		double Npp, double K, double X, int NumGlauberToUse,
		double lowMultCutoff, double highMultCutoff,
		TString outdir, TString OFILE){
//  TFile* glauberFile = new TFile("/Users/matthewharasty/data/Glauber_197_197_32.8mb_WoodsSaxon.root","READ");
//  TFile* multFile = new TFile("/Users/matthewharasty/data/refMult_BES1.root","READ");
//  TH1D* refMult = (TH1D*) multFile->Get("refMult_27");
//  refMult->Rebin(4);

//  GlauberMinimizer* min = new GlauberMinimizer(refMult,glauberFile,100.0);  // last is mult cutoff
//  //min->FitGlauberModel(0.426,1.85,-0.12,100000);
//  min->setChiSqrdMethod("geometric"); // normal / geometric / mixed 
//  min->FitGlauberModel(0.8,2,0.12,-200000); // Npp K X(negative if want fixed at 0.12) NumGlauberToUse(Negative to use All)
	gSystem->Load("../bin/GlauberClass_cxx");
	gSystem->Load("../bin/GlauberMinimizer_cxx");
	
//	cout<<"1"<<endl;
	TFile* glauberFile = new TFile(glauberFileName,"READ");
	TFile* multFile = new TFile(multFileName, "READ");
//	cout<<"2"<<endl;
	TH1D* refMult = (TH1D*) multFile->Get(histname);
	//refMult->Rebin(4)
//	cout<<"3"<<endl;
	GlauberMinimizer* min = new GlauberMinimizer(refMult,glauberFile,lowMultCutoff,highMultCutoff); 
//	cout<<"4"<<endl;
	min->setChiSqrdMethod(chisqmethod); // normal / geometric / mixed 
//	cout<<"5"<<endl;
	min->FitGlauberModel(Npp,K,X,NumGlauberToUse,outdir,OFILE); // Npp K X(negative if want fixed at 0.12) NumGlauberToUse(Negative to use All)
//	cout<<"6"<<endl;
//	return;
}


