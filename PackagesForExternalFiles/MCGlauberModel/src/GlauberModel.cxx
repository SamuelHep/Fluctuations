/* 

   Author: Chris Flores
   Date: July 10, 2014

*/
#include <iostream>
#include <iomanip>

#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TCanvas.h>

#include <GlauberClass.h>

using namespace std;

void GlauberModel(Int_t nEvents, Int_t nNucleonsA, Int_t nNucleonsB,
		  Double_t nnCrossSection, Int_t model=0, TString outDir="data", int seed=-1, Double_t *woodSaxonPars=NULL){

  //OPTIONAL PARAMETERS ARE:
  //model:    changes the model of nucleon distirbution - 0=Hard Sphere, 1=Woods Saxon
  //outDir:   system directory where the output file should be saved. No file will
  //          be saved if outDir=""
  //woodSaxonPars:  An array of wood saxon parameters if you wish to use ones other
  //                than the default values. Useful for systematic error study. Must be of the form
  //                [w,c,a,B20,B40,gamma2,gamma4]
  
  //Create an Instance of our Glauber Class
  GlauberClass glauberEvent(nNucleonsA,nNucleonsB,nnCrossSection);
  
  if ( seed != -1 ) glauberEvent.SetSeed( seed ); 

  //If the user has passed in an empty outDir then no file will be saved otherwise,
  //create the output file name depending on the running configuration
  TString fileName("");
  TFile *outFile = NULL;
  TTree *tree = NULL;
  TBranch *branch = NULL;
  if (outDir.CompareTo("")){
    if ( seed == -1 ) 
      {
	fileName = TString::Format("%s/Glauber_%d_%d_%gmb_%s.root",
				   outDir.Data(),
				   nNucleonsA, nNucleonsB,nnCrossSection,
				   model == 0 ? "UniformSphere" : "WoodsSaxon"
				    );
      }
    else 
      {
	fileName = TString::Format("%s/Glauber_%d_%d_%gmb_%s_%i.root",
				   outDir.Data(),
				   nNucleonsA, nNucleonsB,nnCrossSection,
				   model == 0 ? "UniformSphere" : "WoodsSaxon",
				   seed  );

      }

    outFile = new TFile(fileName,"CREATE");
    tree = new TTree("GlauberTree","GlauberTree");
    branch = tree->Branch("GlauberData",&glauberEvent);
  }
  
  //Loop Over all the Events
  Int_t iEvent(0);
  Int_t XsecTotEvent(0);
  Int_t XsecSucEvent(0);
  int jay=0;
  while (iEvent < nEvents){

    //Set the Max Impact Parameter to 1.3 times larger than the sum
    //of the two nuclear radii
    Double_t maxImpactParameter =
      1.3 * (glauberEvent.GetNuclearRadius(nNucleonsA) +
	     glauberEvent.GetNuclearRadius(nNucleonsB));

    //By default the internally known woods-saxon parameters are used
    //if the user passes in an array of parameters then use it instead
    if (woodSaxonPars)
      glauberEvent.SetWoodsSaxonParameters(woodSaxonPars);

    //Set the impact parameter 
    glauberEvent.SetImpactParameter(0,maxImpactParameter);
    
    //Fill the Nuclei with Nucleons
    glauberEvent.GenerateNucleonPositions(model);
    
    //Count nPart and nColl
    glauberEvent.CountNpartNcoll();
    
    //Fill Tree with events that had at least one participating nucleon
    if (glauberEvent.GetNParticipants() > 0){
      if (tree)
	tree->Fill();
      
      iEvent++;
      XsecSucEvent++;
    }
    XsecTotEvent++;
    
    //Reset the GlauberClass event for reuse
    glauberEvent.Reset();

    //Print the Fraction Complete
    //cout <<setw(15) <<left <<"Percent Complete: " <<(Int_t)(100* (iEvent/(Double_t)nEvents)) <<"\r";
    double perc=100*iEvent/(Double_t)nEvents;
    double percflo=TMath::Floor(perc);
    bool k=false;
    if(jay!=percflo)
      {
        jay=percflo;
        cout<<"Simulation "<<jay<<"% Complete"<<endl;
      }    


  } //End of Loop Over Events
  
  //In case you want to calculate the nuclear cross-section
  Double_t maxImpactParameter =
      1.3 * (glauberEvent.GetNuclearRadius(nNucleonsA) +
             glauberEvent.GetNuclearRadius(nNucleonsB));
  Double_t SuccessRate = 1.0*XsecSucEvent/XsecTotEvent;
  //cout<<"Max Impact Parameter: "<<maxImpactParameter<<endl;
  Double_t TotalArea = TMath::Pi()*(maxImpactParameter)*(maxImpactParameter);
  cout<<"Percent Successful Collisions: "<<100.0*SuccessRate<<endl; 
  cout<<"Cross Section: "<<SuccessRate*TotalArea<<endl;
  //

  if (outFile && branch){
    tree->Write();
    outFile->Close();
  }
  
  cout <<"\nGlauber Model Simulation Complete\n"
       <<"Output file can be found here: " <<fileName <<endl;
}

