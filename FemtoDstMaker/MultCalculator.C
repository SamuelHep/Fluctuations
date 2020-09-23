/*

  MultCalcultor.C
  --Calculates the <Nprotons> at mid rapidity +- 0.5 and generates mult histograms for each centrality bin
  --Use this macro to check the proton distributions after generating the femtoDsts
  
  listname -- pass a list of files as input
  outputname -- name of output file

 */

int Centrality(int fxtmult)
{

  int centBins[8]  = {200,141,118,83,57,38,24,14};

  if (  fxtmult < centBins[0] ) return 0;
  if (  fxtmult < centBins[1] ) return 1;
  if (  fxtmult < centBins[2] ) return 2;
  if (  fxtmult < centBins[3] ) return 3;
  if (  fxtmult < centBins[4] ) return 4;
  if (  fxtmult < centBins[5] ) return 5;
  if (  fxtmult < centBins[6] ) return 6;
  if (  fxtmult < centBins[7] ) return 7;
  return 8;
  
}

void MultCalculator(TString listname = "file.list",TString outputname = "FemtoProtonMult.root")
{

  double mid_rapidity = 1.04;
  TFile * ofile = new TFile(outputname,"recreate");

  gSystem->Load("StFemtoTrack");
  gSystem->Load("StFemtoEvent");

  TH1F * h[9];
  Double_t xbins[8] = {14,24,38,57,83,118,141,200};
  TProfile * mean_prof = new TProfile("mean_profile","mean_profile",7,xbins);

  for( int i=0;i<9;i++)
    {
      h[i] = new TH1F(TString::Format("np_%i",i),"",100,0,100);
    }

  TChain * tc =  new TChain("fDst");

  ifstream filelist(listname);

  string filename;
  while(filelist.good())
    {
      getline(filelist,filename);
      if ( filename.length() != 0 && filename.find(".root") != std::string::npos )
	{
	  tc->Add( filename.c_str() );
	}
    }
  filelist.close();

  StFemtoEvent * event = new StFemtoEvent();
  tc->SetBranchAddress("StFemtoEvent",&event); 

  long int nentries = tc->GetEntries();
  long int tenpercent = nentries/100;
  int percent =0;

  for (int iEntry=0;iEntry<nentries;iEntry++)
    {
      if ( iEntry % tenpercent == 0)
	{
	  cout << percent << "%" << endl;
	  percent += 1;
	}

      tc->GetEntry(iEntry);
      int fxtmult = event->GetFxtMult();
      if (fxtmult < 14) continue;
      
      if ( fabs(event->GetVz() - 200) > 3 ) continue;
      if ( fabs(event->GetVxy()) > 1.5 ) continue;
      
      int Np = 0;

      for (int iTrack=0;iTrack<event->GetEntries();iTrack++)
	{
	  StFemtoTrack trk = event->GetFemtoTrack(iTrack);
	  if ( fabs(trk.GetNSigmaProton()) < 2 )
	    {
	      if ( trk.GetPt() > 0.4 && trk.GetPt() < 2)
		{

		  Double_t mass = 0.938272; // GeV
		  Double_t pt = trk.GetPt();
		  Double_t pz = trk.GetPz();
		  Double_t energy = sqrt( pow(pt,2) + pow(pz,2) + pow(mass,2) );
		  Double_t y = 0.5 * log( ( energy + pz) / (energy - pz) );

		  if ( fabs( fabs( y ) - mid_rapidity ) < 0.5 )
		    {
		      Np += ( trk.GetCharge()>0 ) ? 1 : 0;
		    }
		}
	    }
	}

      mean_prof->Fill(fxtmult,Np);
      h[ Centrality( fxtmult ) ]->Fill( Np );
      
    }


  ofile->cd();
  mean_prof->Write();
  for (int i=0;i<8;i++)
    {
      h[i]->Write();
    }

}
