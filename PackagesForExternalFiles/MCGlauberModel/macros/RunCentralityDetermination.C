#include <iostream>
#include <vector>

using namespace std;

void RunCentralityDetermination(TString DATAFILE, TString DATAHISTO,
				TString GLAUBERFILE,TString OUTFILE,
				Double_t STARTBIN=20, Double_t STOPBIN=-1, 
                                Int_t numTrials=500, Double_t mu1=0.1, Double_t k1=0.1, Double_t hardness=0.12,
                                Double_t npplo=0.00001, Double_t npphi=1.0,
                                Double_t klo=0.00001, Double_t khi=1.0,
                                Double_t xlo=0.0, Double_t xhi=1.0,
                                bool doeverything=1,
                                bool docentralitycuts=1,
				bool freeMu=1, bool freeK=1, bool freeX=1){
 



  cout << "Centrality Determination Settings" << endl;
  cout << "*********************************" << endl;
  cout << "DataFile =" << DATAFILE <<  endl;
  cout << "DataHist =" << DATAHISTO << endl;
  cout << "GlauberFile =" << OUTFILE << endl;
  cout << "OutFile =" << OUTFILE << endl;
  cout << "STARTBin =" << STARTBIN <<  endl;
  cout << "STOPBin =" <<  STOPBIN << endl;
  cout << "number of trials =" <<  numTrials << endl;
  cout << "mu1 = " << mu1 << endl;
  cout << "k1 = " <<  k1 << endl;
  cout << "hardness =" << hardness << endl;
  cout << "nppLow =" << npplo << endl;
  cout << "nppHigh =" <<  npphi << endl;
  cout << "kLow =" <<  klo << endl;
  cout << "kHigh =" << khi << endl;
  cout << "xLo =" << xlo  << endl;
  cout << "xHigh =" << xhi  << endl;
  cout << "Do Everyting =" << doeverything << endl;
  cout << "Do CentCuts = " << docentralitycuts << endl;
  cout << "FreeMu = " << freeMu << endl;
  cout << "FreeK = " << freeK << endl;
  cout << "FreeX = " << freeX << endl;

  //Load the Necessary Libraries
  gSystem->Load("../bin/GlauberClass_cxx");
  gSystem->Load("../bin/GlauberUtil_cxx");
  gSystem->Load("../bin/GlauberModel_cxx");
  gSystem->Load("../bin/FindBestFitNegativeBinomialParameters_cxx");
  gSystem->Load("../bin/FindCentralityBinCuts_cxx");
  gSystem->Load("../bin/DetermineCentralityBin_cxx");
  gSystem->Load("../bin/ReturnStructs_cxx");
  gSystem->Load("../bin/FindNpartNcollDistributions_cxx");
  gSystem->Load("../bin/SystematicErrors_cxx");
  gSystem->Load("../bin/PrintResults_cxx");
  gSystem->Load("../bin/MakeResultsTable_cxx");

  cout<<"Data File: "<<DATAFILE<<endl;
  //Load the Data Histogram that we are trying to match
  TFile *dataFile = new TFile(DATAFILE);
  cout<<"************************************" <<endl;
  if(freeX==0)
  {
  cout<<"           Hardness: "<< hardness <<endl;
  }
  if(freeK==0)
  {
  cout<<"          k Parameter: "<< k1 <<endl;
  }
  if(freeMu==0)
  {
  cout<<"             Npp: "<< mu1 <<endl;
  }
  cout<<"************************************" <<endl;
  TH1D *dataHisto = (TH1D *)dataFile->Get(DATAHISTO);
  TAxis *xaxis = dataHisto->GetXaxis();
  Int_t binnumberstart = xaxis->FindBin(STARTBIN);
  Int_t binnumberstop = xaxis->FindBin(STOPBIN);
  Int_t numberoffitparameters = 3; // 1 for x, 1 for mu, 1 for k
  Int_t numberofdegreesoffreedom = (binnumberstop-binnumberstart+1) - numberoffitparameters; 

  //Load the File from the Glauber Simulation
  GlauberClass *glauberEvent = 0;
  TFile *glauberFile = new TFile(GLAUBERFILE);
  TTree *glauberTree = (TTree *)glauberFile->Get("GlauberTree");
  glauberTree->FindBranch("GlauberData")->SetAddress(&glauberEvent);
  glauberTree->GetEntry(0);
  
  //Create an output file to save things to
  TFile *outFile = new TFile(OUTFILE,"RECREATE");
  
  //Variables Used in the Loops Below
  Double_t startBinCenter(STARTBIN);    //The value of the bin center where we want to start matching
  Double_t stopBinCenter(STOPBIN);      //The value of the bin center where we want to stop matching
  Int_t nTrials(numTrials);                   //The number of choices of (npp,k) permitted
  Double_t Mu(mu1);
  Double_t K(k1);
  Double_t Hardness(hardness);
  Int_t nEvents = glauberTree->GetEntries();
  Bool_t useTrackingEfficiency = true;
  NegBinomialSearchResults bestBinomialParameters;
  Double_t Npplo(npplo);
  Double_t Npphi(npphi);
  Double_t Klo(klo);
  Double_t Khi(khi);
  Double_t Xlo(xlo);
  Double_t Xhi(xhi);

  //Check to make sure we have everything we need
  if (!dataHisto){
    cout <<"ERROR - RunCentralityDetermination.C -- Data histogram not found! ABORTING!" <<endl;
    exit (EXIT_FAILURE);
  }
  if (!glauberTree){
    cout <<"ERROR - RunCentralityDetermination.C -- Glauber Tree not found! ABORTING!" <<endl;
    exit (EXIT_FAILURE);
  }
  
  //________________________________________________________________________________________________
  //First, we need to find the parameters of the negative binomial distribution
  //which results in a multiplicity distrubution that best matches the data.

  cout <<"*** Step 1 of 5: Searching for the Best Negative Binomial Parameters.... " <<endl;
  Int_t numFree = 0;
  if(freeMu==1){numFree++;}
  if(freeK==1){numFree++;}
  if(freeX==1){numFree++;}

  bestBinomialParameters =
    FindBestFitNegativeBinomialParameters(dataHisto,glauberFile,
					  startBinCenter,stopBinCenter,
					  useTrackingEfficiency,nTrials,Mu,K,Hardness,
                                          Npplo, Npphi,
                                          Klo,   Khi,
					  Xlo,   Xhi,
                                          freeMu, freeK, freeX);

  cout <<"  Complete!" <<endl;
  cout <<"  Inverse Chi2: " <<bestBinomialParameters.InverseChi2
       <<" Chi2: " <<1.0/(bestBinomialParameters.InverseChi2)
       <<" Chi2/DOF: "<<(1.0/(bestBinomialParameters.InverseChi2))/(1.0*numberofdegreesoffreedom)
       <<" npp = " <<bestBinomialParameters.npp
       <<" k = "   <<bestBinomialParameters.k 
       <<" x = "   <<bestBinomialParameters.hardness<<endl;
  

  //Save the Default information to the file
  outFile->cd();
  (bestBinomialParameters.bestFitHisto)->SetName("bestFit");
  if(numFree==1){(bestBinomialParameters.chi2Graph1D)->SetName("chiSquaredGraph");}
  if(numFree==2){(bestBinomialParameters.chi2Graph2D)->SetName("chiSquaredGraph");}
  if(numFree==3){(bestBinomialParameters.nTup3D)->SetName("chiSquaredNtup");}
  (bestBinomialParameters.bestFitHisto)->Write();
  if(numFree==1){(bestBinomialParameters.chi2Graph1D)->Write();}
  if(numFree==2){(bestBinomialParameters.chi2Graph2D)->Write();}
  if(numFree==3){(bestBinomialParameters.nTup3D)->Write();}
  dataHisto->Write();

  //Zach added this step to manage crashes
  if(numFree==3)
  {
	TString SUCCESSFILE = "SuccessfulRun_DeleteThisNow.root";
	TFile *successFile = new TFile(SUCCESSFILE,"RECREATE");
  }
  //
  if(docentralitycuts==0 && doeverything==0)
  {
      outFile->cd();
      outFile->Close();
      cout <<"FINISHED SUCCESSFULLY." <<endl;
      return;
  }

  //_______________________________________________________________________________________________
  //Second, Find the centrality bin cuts

  cout <<"*** Step 2 of 5: Determining the centrality bin cuts... " <<endl;
  FindCentralityBinCuts(&bestBinomialParameters);
  cout <<"  Complete!" <<endl;

  for (Int_t iCentBin=0; iCentBin<bestBinomialParameters.centralityBinCuts.size(); iCentBin++){
    cout <<Form("Percent Centrality: %.02F \t Cut: %.05F",
		bestBinomialParameters.centralityBinDefinitions.at(iCentBin),
		bestBinomialParameters.centralityBinCuts.at(iCentBin))
	 <<endl;
  } 

//  if(doeverything==0)
//  {
//      outFile->cd();
//      outFile->Close();
//      cout <<"FINISHED SUCCESSFULLY." <<endl;
//      return;
//  }

  //_______________________________________________________________________________________________
  //Third, Make the Npart, Ncoll, and Impact Parameter Histograms for all the centrality Bins

  cout <<"*** Step 3 of 5: Constructing the Npart, Ncoll, and Impact Parameter Distributions... " <<endl;
  FindNpartNcollDistributions(&bestBinomialParameters,glauberTree,useTrackingEfficiency);
  cout <<"  Complete!" <<endl;
  
  
  if(doeverything==0)
  {
      outFile->cd();
      outFile->mkdir("GlauberHistograms");
      outFile->cd("GlauberHistograms");

  	bestBinomialParameters.nPartTotalHisto->SetName("nPartTotal");
  	bestBinomialParameters.nCollTotalHisto->SetName("nCollTotal");
  	bestBinomialParameters.impactParamTotalHisto->SetName("impactParamTotal");

  	bestBinomialParameters.nPartTotalHisto->Write();
  	bestBinomialParameters.nCollTotalHisto->Write();
  	bestBinomialParameters.impactParamTotalHisto->Write();

      for (Int_t iCentBin=0; iCentBin<bestBinomialParameters.nCentralityBins; iCentBin++){

    	bestBinomialParameters.nPartHistos[iCentBin].SetName(Form("nPart_cent%d",iCentBin));
    	bestBinomialParameters.nCollHistos[iCentBin].SetName(Form("nColl_cent%d",iCentBin));
    	bestBinomialParameters.impactParamHistos[iCentBin].SetName(Form("impactParam_cent%d",iCentBin));

    	bestBinomialParameters.nPartHistos[iCentBin].Write();
    	bestBinomialParameters.nCollHistos[iCentBin].Write();
    	bestBinomialParameters.impactParamHistos[iCentBin].Write();

      }
      outFile->cd();
      outFile->Close();
      cout <<"FINISHED SUCCESSFULLY." <<endl;
      return;
  }
  //_______________________________________________________________________________________________
  //Fourth, Do the systematic Error Investigation

  cout <<"*** Step 4 of 5: Determining the Systematic Errors...." <<endl;
  SystematicErrors(/*nEvents8*/ 1000,glauberEvent->GetNNucleonsNucA(),glauberEvent->GetNNucleonsNucB(),
		   glauberEvent->GetInelasticCrossSection(),&bestBinomialParameters,
		   dataHisto,useTrackingEfficiency);
  cout <<"  Complete!" <<endl;

  //_______________________________________________________________________________________________
  //Fifth, Print Results and Save

  cout <<"*** Step 5 of 5: Printing Results and Saving" <<endl;
  for (Int_t iCentBin=0; iCentBin<bestBinomialParameters.nCentralityBins; iCentBin++){
    PrintResults(&bestBinomialParameters,iCentBin);
  }

  //Make the Results Table
  MakeResultsTable(&bestBinomialParameters);
  
  //Save
  outFile->cd();
  const int nCentBins = bestBinomialParameters.nCentralityBins;
  outFile->WriteObject(&TArrayD(1,&bestBinomialParameters.npp),"npp");
  outFile->WriteObject(&TArrayD(1,&bestBinomialParameters.k),"k");
  outFile->WriteObject(&TArrayD(1,&bestBinomialParameters.InverseChi2),"inverseChi2");
  outFile->WriteObject(&TArrayI(1,&bestBinomialParameters.nCentralityBins),"nCentralityBins");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.centralityBinDefinitions.at(0)),"centralityBinDefinitions");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.centralityBinCuts.at(0)),"centralityBinCuts");

  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nPartMeans.at(0)),"nPartMeans");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nPartStatErrors.at(0)),"nPartStatErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nPartSysErrors.at(0)),"nPartSysErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nCollMeans.at(0)),"nCollMeans");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nCollStatErrors.at(0)),"nCollStatErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.nCollSysErrors.at(0)),"nCollSysErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.impactParamMeans.at(0)),"impactParamMeans");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.impactParamStatErrors.at(0)),"impactParamStatErrors");
  outFile->WriteObject(&TArrayD(nCentBins,&bestBinomialParameters.impactParamSysErrors.at(0)),"impactParamSysErrors");

  outFile->mkdir("GlauberHistograms");
  outFile->cd("GlauberHistograms");

  bestBinomialParameters.nPartTotalHisto->SetName("nPartTotal");
  bestBinomialParameters.nCollTotalHisto->SetName("nCollTotal");
  bestBinomialParameters.impactParamTotalHisto->SetName("impactParamTotal");
  
  bestBinomialParameters.nPartTotalHisto->Write();
  bestBinomialParameters.nCollTotalHisto->Write();
  bestBinomialParameters.impactParamTotalHisto->Write();
  
  for (Int_t iCentBin=0; iCentBin<bestBinomialParameters.nCentralityBins; iCentBin++){

    bestBinomialParameters.nPartHistos[iCentBin].SetName(Form("nPart_cent%d",iCentBin));
    bestBinomialParameters.nCollHistos[iCentBin].SetName(Form("nColl_cent%d",iCentBin));
    bestBinomialParameters.impactParamHistos[iCentBin].SetName(Form("impactParam_cent%d",iCentBin));
    
    bestBinomialParameters.nPartHistos[iCentBin].Write();
    bestBinomialParameters.nCollHistos[iCentBin].Write();
    bestBinomialParameters.impactParamHistos[iCentBin].Write();
    
  }

  outFile->cd();
  outFile->mkdir("Canvases");
  outFile->cd("Canvases");

  bestBinomialParameters.nPartTotalSysErrCanvas->SetName("nPartTotalSysErrCanvas");
  bestBinomialParameters.nCollTotalSysErrCanvas->SetName("nCollTotalSysErrCanvas");
  bestBinomialParameters.impactParamTotalSysErrCanvas->SetName("impactParamTotalSysErrCanvas");
  bestBinomialParameters.resultsTableCanvas->SetName("resultsTableCanvas");
  
  bestBinomialParameters.nPartTotalSysErrCanvas->Write();
  bestBinomialParameters.nCollTotalSysErrCanvas->Write();
  bestBinomialParameters.impactParamTotalSysErrCanvas->Write();
  bestBinomialParameters.resultsTableCanvas->Write();
  
  //Close the File
  outFile->cd();
  outFile->Close();


  cout <<"FINISHED SUCCESSFULLY." <<endl;
  return;
  
}
