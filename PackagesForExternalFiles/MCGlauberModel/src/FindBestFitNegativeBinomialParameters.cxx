#include <iostream>
#include <utility>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TThread.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMutex.h>

#include "GlauberUtil.h"
#include "GlauberClass.h"
#include "ReturnStructs.h"
#include "GlauberUtil.h"

//Globals
TNtuple *nTup3D;
//TGraph2D *chi2Graph2Dhard0;
//TGraph2D *chi2Graph2Dhard1;
//TGraph2D *chi2Graph2Dhard2;
//TGraph2D *chi2Graph2Dhard3;
//TGraph2D *chi2Graph2Dhard4;
//TGraph2D *chi2Graph2Dhard5;
//TGraph2D *chi2Graph2Dhard6;
//TGraph2D *chi2Graph2Dhard7;
//TGraph2D *chi2Graph2Dhard8;
//TGraph2D *chi2Graph2Dhard9;

TH3D *chi3D;
TGraph2D *chi2Graph2D;
TGraph *chi2Graph1D;
TH1D *Prob1D;
Int_t normStartBin;
Int_t normStopBin;
Int_t nGlauberEntries;
Int_t *nPartArr;
Int_t *nCollArr;

TMutex dataHistoMutex;
TMutex nbdHistoMutex;
TMutex chi2Graph1DMutex;
TMutex chi2Graph2DMutex;
TMutex Prob1DMutex;

Bool_t useTrackingEfficiencyGl;


TH1D *MakeMultiplicityForNB(TH1D *dataHisto,TFile *glauberFile,
			    Double_t startMatchingBinCenter,
			    Double_t stopMatchingBinCenter,
			    Double_t npp, Double_t k, Double_t hardness);

struct argStruct{
  TH1D *dataHisto;
  Double_t npp;
  Double_t k;
  Double_t hardness;
  Double_t chi2;
  Int_t index;
  Bool_t freeMu;
  Bool_t freeK;
  Bool_t freeX;
  //Int_t whichhist;
    
};

//________________________________________________________________
//std::pair<double,double> GetNppRange(int numFree){
  
//  const Int_t n 
//  if(numFree==1){n = chi2Graph1D->GetN();}
//  if(numFree==2){n = chi2Graph2D->GetN();}
//  const Int_t top10percent = TMath::CeilNint(n*.10);

//  Int_t index[n];
//  if(numFree==1){TMath::Sort(n,chi2Graph1D->GetY(),index);}
//  if(numFree==2){TMath::Sort(n,chi2Graph2D->GetZ(),index);}

//  Double_t nppOfTops[top10percent];
//  for (Int_t i=0; i<top10percent; i++)
//  {
//    if(numFree==1){nppOfTops[i] = chi2Graph1D->GetX()[index[i]];}
//    if(numFree==2){nppOfTops[i] = chi2Graph2D->GetX()[index[i]];}
//  }

//  Double_t nppMin = TMath::MinElement(top10percent,nppOfTops);
//  Double_t nppMax = TMath::MaxElement(top10percent,nppOfTops);

//  return std::make_pair(nppMin,nppMax);


//}

//________________________________________________________________
//std::pair<double,double> GetKRange(){

//  const Int_t n = chi2Graph->GetN();
//  const Int_t top10percent = TMath::CeilNint(n*.10);

//  Int_t index[n];
//  TMath::Sort(n,chi2Graph->GetZ(),index);

//  Double_t kOfTops[top10percent];
//  for (Int_t i=0; i<top10percent; i++){
//    kOfTops[i] = chi2Graph->GetY()[index[i]];
//  }

//  Double_t kMin = TMath::MinElement(top10percent,kOfTops);
//  Double_t kMax = TMath::MaxElement(top10percent,kOfTops);;

//  return std::make_pair(kMin,kMax);
//}

//_________________________________________________________________
std::vector<double> GetBestPoint(Int_t numfree){
 
  const Int_t n1 = chi2Graph1D->GetN();
  const Int_t n2 = chi2Graph2D->GetN();
  const Int_t n3 = nTup3D->GetEntries();
  vector<double> mkpair;

  
    if(numfree==1)
    {

	Int_t index[n1];
        TMath::Sort(n1,chi2Graph1D->GetY(),index);        
	mkpair.push_back(chi2Graph1D->GetX()[index[0]]);
        mkpair.push_back(chi2Graph1D->GetY()[index[0]]);
    }
    if(numfree==2)
    {
	Int_t index[n2];
        TMath::Sort(n2,chi2Graph2D->GetZ(),index);
        mkpair.push_back(chi2Graph2D->GetX()[index[0]]);
        mkpair.push_back(chi2Graph2D->GetY()[index[0]]);
        mkpair.push_back(chi2Graph2D->GetZ()[index[0]]);
        //cout<<"List of chis: "<<chi2Graph2D->GetZ()[index[0]]<<chi2Graph2D->GetZ()[index[1]]<<chi2Graph2D->GetZ()[index[2]]<<chi2Graph2D->GetZ()[index[3]]<<chi2Graph2D->GetZ()[index[4]]<<endl;
        //cout<<"BestPintki: "<<mkpair[2]<<endl;
    }
    if(numfree==3)
    {
        Int_t index[n3];   
	//TMath::Sort(n3,nTup3D->GetA(),index); //Fix This
        //mkpair.push_back(nTup3D->GetX()[index[0]]);
        //mkpair.push_back(nTup3D->GetY()[index[0]]);
        //mkpair.push_back(nTup3D->GetZ()[index[0]]);
        //mkpair.push_back(nTup3D->GetA()[index[0]]);
    

        Float_t Npp;
        Float_t K;
        Float_t Hardness;
        Float_t InverseChiSquared;
        nTup3D->SetBranchAddress("Npp",&Npp);
        nTup3D->SetBranchAddress("K",&K);
        nTup3D->SetBranchAddress("Hardness",&Hardness);
        nTup3D->SetBranchAddress("InverseChiSquared",&InverseChiSquared);
        
        Double_t thismu=0.0;
        Double_t thisk=0.0;
        Double_t thishardness=0.0;
        Double_t thischi2=0.0;

        for(int entry=0;entry<n3;entry++)
        {
            nTup3D->GetEvent(entry);
            if(InverseChiSquared>thischi2)
	    {
		thismu=Npp;
		thisk=K;
		thishardness=Hardness;
		thischi2=InverseChiSquared;
	    }
        }
        mkpair.push_back(thismu);
	mkpair.push_back(thisk);
	mkpair.push_back(thishardness);
	mkpair.push_back(thischi2);
     }

  return mkpair;
}

//_________________________________________________________________
Double_t ComputeHistogramChiSquared(TH1D *hData, TH1D *hSim, Int_t startBin, Int_t stopBin){

  Double_t chi2(0);

  for (Int_t iBin=startBin; iBin<stopBin; iBin++){

    //Don't consider points from data histo that have no error
    if (hData->GetBinError(iBin) == 0.0)
      continue;
    
    // Uncomment the following for log chi-squared
    //
   // if(hSim->GetBinContent(iBin)==0)
   // {
  //  	chi2 += 100;
   // }
    //
    //Double_t num = 50.0*(TMath::Abs(TMath::Log10(hData->GetBinContent(iBin))) - TMath::Abs(TMath::Log10(hSim->GetBinContent(iBin))));
    //Double_t denom = 1.0;
    //Double_t denom = TMath::Abs(TMath::Log10(hData->GetBinError(iBin)));
    //
    // Uncomment the following for linear chi-squared
    //
    Double_t num = (hData->GetBinContent(iBin)) - (hSim->GetBinContent(iBin));    
    Double_t denom = (hData->GetBinError(iBin));
   
    //chi2 += iBin*pow(num/denom,2);
    chi2 += pow(num/denom,2);
  }
  if((1.0/chi2)==0.0)
 {
     return 1000000000.0;    
  }
  //if(1.0/chi2==0.0){cout<<"invchi2: "<<1.0/chi2<<endl;}
  return chi2;
}

//___________________________________________________________________
void *DoParticleProduction(void *args){

  argStruct *argVals = (argStruct *)args;
  Double_t npp = argVals->npp;
  Double_t k   = argVals->k;
  Double_t hardness = argVals->hardness;
  Int_t index  = argVals->index;
  Bool_t freeMu = argVals->freeMu;
  Bool_t freeK = argVals->freeK;
  Bool_t freeX = argVals->freeX;
  //Int_t whichhist = argVals->whichhist;

  //Create an instance of the GlauberClass
  GlauberClass glauberEvent;
  glauberEvent.SetNegativeBinomialParameters(npp,k);
  nbdHistoMutex.Lock();
  TH1D * nbdHist = MakeNegativeBinomialHist(npp,k,hardness,index,glauberEvent);
  nbdHistoMutex.UnLock();

  dataHistoMutex.Lock();
  TH1D dataHisto(*argVals->dataHisto);
  dataHistoMutex.UnLock();

  Double_t nParticles(0);
  TH1D *simHisto = new TH1D(Form("hTemp_%d",index),Form("hTemp_%d",index),
			 dataHisto.GetNbinsX(),
			 dataHisto.GetBinLowEdge(1),
			 dataHisto.GetBinLowEdge(dataHisto.GetNbinsX())+dataHisto.GetBinWidth(1));

  for (Int_t i=0; i<nGlauberEntries; i++){
    nParticles = glauberEvent.ProduceParticles(nPartArr[i],nCollArr[i],nbdHist,hardness,useTrackingEfficiencyGl);
    simHisto->Fill(nParticles);
//    if(nParticles > 0)
//	{
//    	AVERAGENPPVECTOR.push_back(nParticles);
//        }
  }

  delete nbdHist;

  //Prepare the simulated histo for comparison to data histo
  simHisto->Sumw2();
  Double_t scaleFactor = simHisto->Integral(normStartBin,simHisto->GetNbinsX());

  //If this simHisto has a bad scale factor then just return
  if (scaleFactor < 1){
    TThread::Lock();
    cout <<"Thread " <<index <<" had a bad scale factor" <<endl;
    TThread::UnLock();
    return (void *)NULL;
  }

  //Scale the Histogram
  simHisto->Scale(1.0/(Double_t)simHisto->Integral(normStartBin,simHisto->GetNbinsX()));

//  Double_t AVERAGENPPSIZE = AVERAGENPPVECTOR.size();
//  Double_t AVERAGENPPSUM = 0;
//  for (int num=0; num<AVERAGENPPSIZE; num++)
//	{
//	AVERAGENPPSUM += 
//	}
//  Double_t AVERAGENPP;
//
 
  
  Double_t chi2 = ComputeHistogramChiSquared(&dataHisto,simHisto,normStartBin,normStopBin);
  //if(1.0/chi2==0.0){cout<<"invchiinpartcileproduction"<<1.0/chi2<<endl;}
  chi2Graph1DMutex.Lock();
 // chi2Graph2DMutex.Lock();
 // Prob1DMutex.Lock();
    if(freeMu==1 && freeK==1 && freeX==0){chi2Graph2D->SetPoint(chi2Graph2D->GetN(),npp,k,1.0/chi2);}
    if(freeMu==0 && freeK==1 && freeX==1){chi2Graph2D->SetPoint(chi2Graph2D->GetN(),k,hardness,1.0/chi2);}
    if(freeMu==1 && freeK==0 && freeX==1){chi2Graph2D->SetPoint(chi2Graph2D->GetN(),npp,hardness,1.0/chi2);}
    if(freeMu==1 && freeK==0 && freeX==0){chi2Graph1D->SetPoint(chi2Graph1D->GetN(),npp,1.0/chi2);Prob1D->Fill(npp,1.0/chi2);}
    if(freeMu==0 && freeK==1 && freeX==0){chi2Graph1D->SetPoint(chi2Graph1D->GetN(),k,1.0/chi2);Prob1D->Fill(k,1.0/chi2);}
    if(freeMu==0 && freeK==0 && freeX==1){chi2Graph1D->SetPoint(chi2Graph1D->GetN(),hardness,1.0/chi2);Prob1D->Fill(hardness,1.0/chi2);}
    if(freeMu==1 && freeK==1 && freeX==1)
    { 
          // if(whichhist==0){chi2Graph2Dhard0->SetPoint(chi2Graph2Dhard0->GetN(),npp,k,1.0/chi2);}
          // else if(whichhist==1){chi2Graph2Dhard1->SetPoint(chi2Graph2Dhard1->GetN(),npp,k,1.0/chi2);}
          // else if(whichhist==2){chi2Graph2Dhard2->SetPoint(chi2Graph2Dhard2->GetN(),npp,k,1.0/chi2);}
          // else if(whichhist==3){chi2Graph2Dhard3->SetPoint(chi2Graph2Dhard3->GetN(),npp,k,1.0/chi2);}
          // else if(whichhist==4){chi2Graph2Dhard4->SetPoint(chi2Graph2Dhard4->GetN(),npp,k,1.0/chi2);}
          // else if(whichhist==5){chi2Graph2Dhard5->SetPoint(chi2Graph2Dhard5->GetN(),npp,k,1.0/chi2);}
          // else if(whichhist==6){chi2Graph2Dhard6->SetPoint(chi2Graph2Dhard6->GetN(),npp,k,1.0/chi2);}
          // else if(whichhist==7){chi2Graph2Dhard7->SetPoint(chi2Graph2Dhard7->GetN(),npp,k,1.0/chi2);}           
          // else if(whichhist==8){chi2Graph2Dhard8->SetPoint(chi2Graph2Dhard8->GetN(),npp,k,1.0/chi2);}
          // else if(whichhist==9){chi2Graph2Dhard9->SetPoint(chi2Graph2Dhard9->GetN(),npp,k,1.0/chi2);}
           
           chi3D->Fill(npp,k,hardness,1.0/chi2);
           nTup3D->Fill(npp,k,hardness,1.0/chi2);
    }
 chi2Graph1DMutex.UnLock();
 // chi2Graph2DMutex.UnLock();
 // Prob1DMutex.UnLock();

  delete simHisto;

  return NULL;
}



//______________________________________________________________________________________________
NegBinomialSearchResults
FindBestFitNegativeBinomialParameters(TH1D *dataHisto,TFile *glauberFile,
				      Double_t startMatchingBinCenter,
				      Double_t stopMatchingBinCenter,
				      Bool_t useTrackingEfficiency,
				      const Int_t nThreads=1000,
                                      Double_t npp=0.1, Double_t k=0.1, Double_t hardness=0.5,
                                      Double_t npplo=0.00001, Double_t npphi=1.0,
                                      Double_t klo=0.00001,   Double_t khi=0.1,
                                      Double_t xlo=0.0,       Double_t xhi=1.0,
                                      Bool_t freeMu=true, Bool_t freeK=true, Bool_t freeX=true){

  useTrackingEfficiencyGl = useTrackingEfficiency;
   
  Double_t nppMin(npplo); //0.1
  Double_t nppMax(npphi); //2.0
  Double_t kMin(klo);   //0.0
  Double_t kMax(khi);     //5
  Double_t xMin(xlo);
  Double_t xMax(xhi);
  
  
  //The bin for which will will start the matchin procedure
  normStartBin = dataHisto->FindBin(startMatchingBinCenter);

  //This is the bin for which we will end the matching procedure.
  //By default the last bin of the dataHisto is used as the stopping bin.
  //However, if the user has passed in a valid value for the stop bin then use it instead.
  normStopBin = dataHisto->GetNbinsX();
  if (stopMatchingBinCenter > startMatchingBinCenter)
    normStopBin = dataHisto->FindBin(stopMatchingBinCenter);

  //Build the Npart and Ncoll Arrays from the file so that we can do the
  //particle production procedure using multiple threads
  GlauberClass *glauberEvent = 0;
  TTree *glauberTree = (TTree *)glauberFile->Get("GlauberTree");
  glauberTree->FindBranch("GlauberData")->SetAddress(&glauberEvent);
  
  const Int_t nEntries = glauberTree->GetEntries();
  nGlauberEntries = nEntries;
  nPartArr = new Int_t[nEntries];
  nCollArr = new Int_t[nEntries];

  for (Int_t i=0; i<glauberTree->GetEntries(); i++){
    glauberTree->GetEntry(i);
    nPartArr[i] = glauberEvent->GetNParticipants();
    nCollArr[i] = glauberEvent->GetNBinaryCollisions();
  }//End Loop Over Glauber Tree


  //Make sure we have everything we need
  if (!dataHisto){
    cout <<"ERROR: FindBestFitNegativeBinomialParameters() - dataHisto does not exist!" <<endl;
    exit (EXIT_FAILURE);
  }

  //Prepare the Data Histogram...do this only if sumw2 hasn't been called before
  //this is to protect the histogram in the case of mulitple iterations
  if ( ((TArrayD *)dataHisto->GetSumw2())->GetSize() == 0){
    dataHisto->Sumw2();
  }
  dataHisto->Scale(1.0/(Double_t)dataHisto->Integral(normStartBin,normStopBin));
  dataHisto->SetMarkerColor(kBlack);
  dataHisto->SetMarkerStyle(kFullCircle);
  
  TRandom3 rand(0);
    
    // This is where the new types of scans are being implemented
    Int_t numFree=0;
    if(freeMu==1){numFree++;}
    if(freeK==1){numFree++;}
    if(freeX==1){numFree++;}
    
    //if(numFree==3)
    //{
    //    cout <<"ERROR: You're asking me to do too much. I can't do a three dimensional chi-squared minimization, so either fix one of the free parameters or rewrite the code." <<endl;
    //    exit (EXIT_FAILURE);
    //}
    if(numFree==0)
    {
        cout <<"ERROR: You're asking me to do too little. I can't do a 0 dimensional chi-squared minimization, so either free one or two of the parameters or rewrite the code." <<endl;
        exit (EXIT_FAILURE);
    }
    nTup3D = new TNtuple("ntup","Glauber Model Data","Npp:K:Hardness:InverseChiSquared");
    
    chi3D = new TH3D("chi3D","Chi Squared Graph", 20, nppMin, nppMax, 20, kMin, kMax, 20, xMin, xMax); chi3D->SetMarkerStyle(kDot);
    //chi2Graph2Dhard0 = new TGraph2D(); chi2Graph2Dhard0->SetMarkerStyle(kDot);
    //chi2Graph2Dhard1 = new TGraph2D(); chi2Graph2Dhard1->SetMarkerStyle(kDot);
    //chi2Graph2Dhard2 = new TGraph2D(); chi2Graph2Dhard2->SetMarkerStyle(kDot);
    //chi2Graph2Dhard3 = new TGraph2D(); chi2Graph2Dhard3->SetMarkerStyle(kDot);
    //chi2Graph2Dhard4 = new TGraph2D(); chi2Graph2Dhard4->SetMarkerStyle(kDot);
    //chi2Graph2Dhard5 = new TGraph2D(); chi2Graph2Dhard5->SetMarkerStyle(kDot);
    //chi2Graph2Dhard6 = new TGraph2D(); chi2Graph2Dhard6->SetMarkerStyle(kDot);
    //chi2Graph2Dhard7 = new TGraph2D(); chi2Graph2Dhard7->SetMarkerStyle(kDot);
    //chi2Graph2Dhard8 = new TGraph2D(); chi2Graph2Dhard8->SetMarkerStyle(kDot);
    //chi2Graph2Dhard9 = new TGraph2D(); chi2Graph2Dhard9->SetMarkerStyle(kDot);
   
   
    chi2Graph2D = new TGraph2D(); chi2Graph2D->SetMarkerStyle(kDot);
    chi2Graph1D = new TGraph();   chi2Graph1D->SetMarkerStyle(kFullCircle); 
   if(numFree==1)
	{
	if(freeMu==1){Prob1D = new TH1D("Prob1D", "Probability Dist", 20, nppMin, nppMax);}
	if(freeK==1){Prob1D = new TH1D("Prob1D", "Probability Dist", 20, kMin, kMax);}
	if(freeX==1){Prob1D = new TH1D("Prob1D", "Probability Dist", 20, xMin, xMax);}
	}
		
  
  const Int_t nLoops(3);
  const Int_t nConcurrentThreads(25);
  double jay=0;
  int eye=0;
  for (Int_t iLoop=0; iLoop<nLoops; iLoop++){

    Int_t nSubmittedThreads(0);

    //Decrease the number of threads for each loop
    Int_t nThreadsToSubmit = nThreads/(double)(iLoop+1);

    //While the number of running threads is less than the
    //number of number of threads to be submitted for this loop...
    while (nSubmittedThreads < nThreadsToSubmit){
      double perc=100.0*eye/(nThreads/(double)(1)+nThreads/(double)(2)+nThreads/(double)(3));
      double percflo=TMath::Floor(perc);
      if(jay!=percflo)
      {
        jay=percflo;
        cout<<"Calculation "<<jay<<"% Complete"<<endl;
      }

      //Only submit more threads when there are less than
      //nConcurrentThreads running
      std::vector <argStruct *> argVector;
      std::vector <TThread *>   threadVector;
      while (TThread::Exists() < nConcurrentThreads){

	argVector.push_back(new argStruct);
	argVector.back()->index = nSubmittedThreads;
	argVector.back()->dataHisto = dataHisto;
    
    if(freeMu==0){argVector.back()->npp = npp;}
    if(freeK==0){argVector.back()->k = k;}
    if(freeX==0){argVector.back()->hardness = hardness;}

	//Set the nb parameters based on which loop we are in
	if (iLoop == 0)
    
	{
        if(freeMu==1){argVector.back()->npp = rand.Uniform(nppMin,nppMax);}
        if(freeK==1){argVector.back()->k   = rand.Uniform(kMin,kMax);}
        if(freeX==1){argVector.back()->hardness   = rand.Uniform(xMin,xMax);}
        //cout<<"npp: "<<rand.Uniform(nppMin,nppMax)<<endl;
	}
	else if (iLoop == 1)
        {   
	  //std::pair<double,double> nppRange = GetNppRange();
        if(freeMu==1){argVector.back()->npp = rand.Uniform(nppMin,nppMax);}
        if(freeK==1){argVector.back()->k   = rand.Uniform(kMin,kMax);}
        if(freeX==1){argVector.back()->hardness   = rand.Uniform(xMin,xMax);}
	//cout<<"npp: "<<rand.Uniform(nppMin,nppMax)<<endl;
	}
       // Int_t whichhist=1000;
       // if(numFree==3)
       // {
       //    whichhist=100;
       //    if(xMin+0.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+1.0*(xMax-xMin)/10.0){whichhist=0;}
       //    else if(xMin+1.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+2.0*(xMax-xMin)/10.0){whichhist=1;}
       //    else if(xMin+2.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+3.0*(xMax-xMin)/10.0){whichhist=2;}
       //    else if(xMin+3.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+4.0*(xMax-xMin)/10.0){whichhist=3;}
       //    else if(xMin+4.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+5.0*(xMax-xMin)/10.0){whichhist=4;}
       //    else if(xMin+5.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+6.0*(xMax-xMin)/10.0){whichhist=5;}
       //    else if(xMin+6.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+7.0*(xMax-xMin)/10.0){whichhist=6;}
       //    else if(xMin+7.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+8.0*(xMax-xMin)/10.0){whichhist=7;}
       //    else if(xMin+8.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+9.0*(xMax-xMin)/10.0){whichhist=8;}
       //    else if(xMin+9.0*(xMax-xMin)/10.0 <= hardness && hardness <= xMin+10.0*(xMax-xMin)/10.0){whichhist=9;}
       //    else if(whichhist==100){cout<<"ERROR!: Splitting the parameter space into slices in hardness didn't work!  >:("<<endl;}
       // }

	else{
	  Double_t tempNpp(0), tempK(0), tempX(0);
	 // Double_t tempNpp, tempK, tempX;
        if(numFree==3)
	{
	//	argVector.back()->hardness   = rand.Uniform(xMin,xMax);
        //        whichhist=100;
        //        if(xMin+0.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+1.0*(xMax-xMin)/10.0){whichhist=0; chi2Graph2Dhard0->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(xMin+1.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+2.0*(xMax-xMin)/10.0){whichhist=1; chi2Graph2Dhard1->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(xMin+2.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+3.0*(xMax-xMin)/10.0){whichhist=2; chi2Graph2Dhard2->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(xMin+3.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+4.0*(xMax-xMin)/10.0){whichhist=3; chi2Graph2Dhard3->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(xMin+4.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+5.0*(xMax-xMin)/10.0){whichhist=4; chi2Graph2Dhard4->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(xMin+5.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+6.0*(xMax-xMin)/10.0){whichhist=5; chi2Graph2Dhard5->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(xMin+6.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+7.0*(xMax-xMin)/10.0){whichhist=6; chi2Graph2Dhard6->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(xMin+7.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+8.0*(xMax-xMin)/10.0){whichhist=7; chi2Graph2Dhard7->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(xMin+8.0*(xMax-xMin)/10.0 <= hardness && hardness < xMin+9.0*(xMax-xMin)/10.0){whichhist=8; chi2Graph2Dhard8->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(xMin+9.0*(xMax-xMin)/10.0 <= hardness && hardness <= xMin+10.0*(xMax-xMin)/10.0){whichhist=9; chi2Graph2Dhard9->GetHistogram()->GetRandom2(tempNpp, tempK);}
        //        else if(whichhist==100){cout<<"ERROR!: Splitting the parameter space into slices in hardness didn't work!  >:("<<endl;}                

                //argVector.back()->npp = rand.Uniform(nppMin,nppMax);
		//argVector.back()->k   = rand.Uniform(kMin,kMax);
	    chi3D->GetRandom3(tempNpp, tempK, tempX);
            argVector.back()->npp = tempNpp;
            argVector.back()->k   = tempK;
            argVector.back()->hardness = tempX;
	
	}
	if(freeMu==1 && freeK==1 && freeX==0)
        {
            chi2Graph2D->GetHistogram()->GetRandom2(tempNpp, tempK);
            argVector.back()->npp = tempNpp;
            argVector.back()->k   = tempK;
        }
        if(freeMu==1 && freeK==0 && freeX==1)
        {
            chi2Graph2D->GetHistogram()->GetRandom2(tempNpp, tempX);
            argVector.back()->npp = tempNpp;
            argVector.back()->hardness   = tempX;
        }
        if(freeMu==0 && freeK==1 && freeX==1)
        {
            chi2Graph2D->GetHistogram()->GetRandom2(tempK, tempX);
            argVector.back()->k = tempK;
            argVector.back()->hardness   = tempX;
        }
        if(freeK==1 && freeX==0 && freeMu==0)
        {   
	    //chi2Graph1D->GetHistogram()->GetRandom();
            tempK=Prob1D->GetRandom();
	    argVector.back()->k = tempK;
        }
        if(freeK==0 && freeX==1 && freeMu==0)
        {
            //tempX = chi2Graph1D->GetHistogram()->GetRandom();
            tempX = Prob1D->GetRandom();
	    argVector.back()->hardness = tempX;
        }
        if(freeK==0 && freeX==0 && freeMu==1)
        {
            tempNpp = Prob1D->GetRandom();
            argVector.back()->npp = tempNpp;
       	 //   cout<<"npp: "<<tempNpp<<endl;
	}
	}

    argVector.back()->freeMu = freeMu;
    argVector.back()->freeK = freeK;
    argVector.back()->freeX = freeX;
    //argVector.back()->whichhist = whichhist;
	//Create the new thread and then run it
	threadVector.push_back(new TThread(Form("thread_%d",nSubmittedThreads),DoParticleProduction,
					   (void *)argVector.back()));
	threadVector.back()->Run();

	//Incremenet the number of threads that have been submitted
      eye++;
      nSubmittedThreads++;


      }//End Loop Over submitting nConcurrent Threads

      //After submitting the nConcurrent Threads wait for them to be done
      //before starting the next batch
      for (unsigned int i=0; i<threadVector.size(); i++){
	threadVector.at(i)->Join();
	threadVector.at(i)->Delete();
	delete argVector.at(i);
      }

    }//End Loop Over submitted Threads

  }//End Loops

  //Best Npp and K and X
  
  std::vector<double> bestPoint = GetBestPoint(numFree);
  //cout<<"BestPoint2: "<<bestPoint[3];
  NegBinomialSearchResults results;
    
    if(freeMu==1 && freeK==1 && freeX==1)
    {
        results.npp = bestPoint[0];
        results.k = bestPoint[1];
        results.hardness = bestPoint[2];
        results.InverseChi2 = bestPoint[3];
        results.bestFitHisto = MakeMultiplicityForNB(dataHisto,glauberFile,
                                                     startMatchingBinCenter,stopMatchingBinCenter,
                                                     results.npp,results.k,results.hardness);
        results.nTup3D = (TNtuple *)nTup3D->Clone(nTup3D->GetName());
    } 
    if(freeMu==1 && freeK==1 && freeX==0)
    {
        results.npp = bestPoint[0];
        results.k = bestPoint[1];
        results.hardness = hardness;
        results.InverseChi2 = chi2Graph2D->Interpolate(results.npp,results.k);
        if(results.InverseChi2==0.0){results.InverseChi2 = bestPoint[2];}
	//cout<<"InverseChi2: "<<results.InverseChi2<<endl;
        results.bestFitHisto = MakeMultiplicityForNB(dataHisto,glauberFile,
                                                     startMatchingBinCenter,stopMatchingBinCenter,
                                                     results.npp,results.k,hardness);
        results.chi2Graph2D = (TGraph2D *)chi2Graph2D->Clone(chi2Graph2D->GetName());
    }
    if(freeMu==0 && freeK==1 && freeX==1)
    {
        results.npp = npp;
        results.k = bestPoint[0];
        results.hardness = bestPoint[1];
        results.InverseChi2 = chi2Graph2D->Interpolate(results.k,results.hardness);
        if(results.InverseChi2==0.0){results.InverseChi2 = bestPoint[2];}
        results.bestFitHisto = MakeMultiplicityForNB(dataHisto,glauberFile,
                                                     startMatchingBinCenter,stopMatchingBinCenter,
                                                     npp,results.k,results.hardness);
        results.chi2Graph2D = (TGraph2D *)chi2Graph2D->Clone(chi2Graph2D->GetName());
    }
    if(freeMu==1 && freeK==0 && freeX==1)
    {
        results.npp = bestPoint[0];
        results.k = k;
        results.hardness = bestPoint[1];
        results.InverseChi2 = chi2Graph2D->Interpolate(results.npp,results.hardness);
        if(results.InverseChi2==0.0){results.InverseChi2 = bestPoint[2];}
        results.bestFitHisto = MakeMultiplicityForNB(dataHisto,glauberFile,
                                                     startMatchingBinCenter,stopMatchingBinCenter,
                                                     results.npp,k,results.hardness);
        results.chi2Graph2D = (TGraph2D *)chi2Graph2D->Clone(chi2Graph2D->GetName());
    }
    if(freeMu==1 && freeK==0 && freeX==0)
    {
        results.npp = bestPoint[0];
        results.k = k;
        results.hardness = hardness;
        results.InverseChi2 = bestPoint[1];
        results.bestFitHisto = MakeMultiplicityForNB(dataHisto,glauberFile,
                                                     startMatchingBinCenter,stopMatchingBinCenter,
                                                     results.npp,k,hardness);
        results.chi2Graph1D = (TGraph *)chi2Graph1D->Clone(chi2Graph1D->GetName());
    }
    if(freeMu==0 && freeK==1 && freeX==0)
    {   
        results.npp = npp;
        results.k = bestPoint[0];
        results.hardness = hardness;
	
        results.InverseChi2 = bestPoint[1];
        results.bestFitHisto = MakeMultiplicityForNB(dataHisto,glauberFile,
                                                     startMatchingBinCenter,stopMatchingBinCenter,
                                                     npp,results.k,hardness);
        
	results.chi2Graph1D = (TGraph *)chi2Graph1D->Clone(chi2Graph1D->GetName());
    	
    }
    if(freeMu==0 && freeK==0 && freeX==1)
    {
        results.npp = npp;
        results.k = k;
        results.hardness = bestPoint[0];
        results.InverseChi2 = bestPoint[1];
        results.bestFitHisto = MakeMultiplicityForNB(dataHisto,glauberFile,
                                                     startMatchingBinCenter,stopMatchingBinCenter,
                                                     npp,k,results.hardness);
        results.chi2Graph1D = (TGraph *)chi2Graph1D->Clone(chi2Graph1D->GetName());
    }
    
 // results.InverseChi2 = bestPoint[2];
 
    //results.InverseChi2 = chi2Graph->Interpolate(results.npp,results.k);
// ////////////////////////////////////////////
// ////////////////////////////////////////////
// ////////////////////////////////////////////
//   const Int_t n = chi2Graph->GetZ(),index);
//   results.InverseChi2 = chi2Graph->GetZ()[index[0]];
//   cout<<"INVERSECHI2: "<<chi2Graph->GetZ()[index[0]]<<endl;
//   cout<<"REAL MIN: "<<chi2Graph->GetZmax()<<endl;
//  results.InverseChi2 = chi2Graph->GetMinimum();
//   for(int aye=0;aye<n;aye++){cout<<"ChiSq: "<<chi2Graph->GetZ()[index[aye]]<<endl;}
// ////////////////////////////////////////////
// ////////////////////////////////////////////
// //////////////////////////////////////////// 
 
  //results.bestFitHisto = MakeMultiplicityForNB(dataHisto,glauberFile,
  //					       startMatchingBinCenter,stopMatchingBinCenter,
  //					       results.npp,results.k,hardness);
  //results.chi2Graph = (TGraph2D *)chi2Graph->Clone(chi2Graph->GetName());

  //Clean Up
  delete nPartArr;
  delete nCollArr;
  delete chi2Graph2D;
  delete chi3D;
  //delete chi2Graph2Dhard0;
  //delete chi2Graph2Dhard1;
  //delete chi2Graph2Dhard2;
  //delete chi2Graph2Dhard3;
  //delete chi2Graph2Dhard3;
  //delete chi2Graph2Dhard4;
  //delete chi2Graph2Dhard5;
  //delete chi2Graph2Dhard6;
  //delete chi2Graph2Dhard7;
  //delete chi2Graph2Dhard8;
  //delete chi2Graph2Dhard9;
  delete chi2Graph1D;
  delete nTup3D;
  delete Prob1D;
  
  return results;

}


//_______________________________________________________________________________________
TH1D *MakeMultiplicityForNB(TH1D *dataHisto,TFile *glauberFile,
			   Double_t startMatchingBinCenter,
			   Double_t stopMatchingBinCenter,
			   Double_t npp, Double_t k, Double_t hardness){

  //The bin for which will will start the matchin procedure
  normStartBin = dataHisto->FindBin(startMatchingBinCenter);

  //This is the bin for which we will end the matching procedure.
  //By default the last bin of the dataHisto is used as the stopping bin.
  //However, if the user has passed in a valid value for the stop bin then use it instead.
  normStopBin = dataHisto->GetNbinsX();
  if (stopMatchingBinCenter > startMatchingBinCenter)
    normStopBin = dataHisto->FindBin(stopMatchingBinCenter);

  //Prepare the Data Histogram...do this only if sumw2 hasn't been called before
  //this is to protect the histogram in the case of mulitple iterations
  if ( ((TArrayD *)dataHisto->GetSumw2())->GetSize() == 0){
    dataHisto->Sumw2();
    dataHisto->Scale(1.0/(Double_t)dataHisto->Integral(normStartBin,dataHisto->GetNbinsX()));
    dataHisto->SetMarkerColor(kBlack);
   // dataHisto->SetMarkerStyle(kFullCircle);
  }
  
  //Build the Npart and Ncoll Arrays from the file so that we can do the
  //particle production procedure using multiple threads
  GlauberClass *glauberEvent = 0;
  TTree *glauberTree = (TTree *)glauberFile->Get("GlauberTree");
  glauberTree->FindBranch("GlauberData")->SetAddress(&glauberEvent);
  glauberEvent->SetNegativeBinomialParameters(npp,k);
  TH1D * nbdHist = MakeNegativeBinomialHist(npp,k,hardness,0,*glauberEvent);

  Double_t nParticles(0);
  TH1D *simHisto = new TH1D("htemp","hTemp",
			    dataHisto->GetNbinsX(),
			    dataHisto->GetBinLowEdge(1),
			    dataHisto->GetBinLowEdge(dataHisto->GetNbinsX())+dataHisto->GetBinWidth(1));

  const Int_t nEntries = glauberTree->GetEntries();
  for (Int_t i=0; i<nEntries; i++){
    glauberTree->GetEntry(i);
    nParticles = glauberEvent->ProduceParticles(hardness,nbdHist, useTrackingEfficiencyGl);
    simHisto->Fill(nParticles);

  }//End Loop Over Glauber Tree

  //Prepare the simulated histo for comparison to data histo
  simHisto->Sumw2();
  Double_t scaleFactor = simHisto->Integral(normStartBin,simHisto->GetNbinsX());
  simHisto->Scale(1.0/scaleFactor);
  simHisto->SetLineColor(kRed);

  return simHisto;
}
