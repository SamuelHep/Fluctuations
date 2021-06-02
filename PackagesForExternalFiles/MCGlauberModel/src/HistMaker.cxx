#include "TTree.h"
#include "TFile.h"

#include "HistMaker.h"
#include "GlauberClass.h"
#include "GlauberUtil.h"


void HistMaker(TString DATAFILE,
	       TString OUTFILE,
	       //	       std::vector<int> centVals,
	       double npp,
	       double k,
	       double hardness
	       )
{
  
  //Load the File from the Glauber Simulation
  GlauberClass *glauberEvent = 0;
  TFile *glauberFile = new TFile(DATAFILE);
  TTree *glauberTree = (TTree *)glauberFile->Get("GlauberTree");
  glauberTree->FindBranch("GlauberData")->SetAddress(&glauberEvent);
  glauberTree->GetEntry(0);

  glauberEvent->SetNegativeBinomialParameters(npp,k);
  //Make NBD
  TH1D * nbdHist = MakeNegativeBinomialHist(npp,k,hardness,0,*glauberEvent);  

  //Create an output file to save things to
  TFile *outFile = new TFile(OUTFILE,"RECREATE");
 
  int nColl=0;
  int nPart=0;
  int nParticles=0;
  double b=0;

  //SimpleTree to fill
  TTree * simpleTree = new TTree("reducedGlauber","");
  simpleTree->Branch("b",&b);
  simpleTree->Branch("nColl",&nColl);
  simpleTree->Branch("nPart",&nPart);
  simpleTree->Branch("nParticles",&nParticles);

  for ( Int_t i=0; i < glauberTree->GetEntries(); i++)
    {
      glauberTree->GetEntry(i);
      b          = glauberEvent->GetImpactParameter();
      nColl      = glauberEvent->GetNBinaryCollisions();
      nPart      = glauberEvent->GetNParticipants();
      nParticles = glauberEvent->ProduceParticles(nPart,nColl,nbdHist,hardness, true);
      simpleTree->Fill();
    }
 
  outFile->cd();
  simpleTree->Write();

}
