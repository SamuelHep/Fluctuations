//Class for proton efficiency
//Currently only has constant efficiency
//should be updated to have efficiency curves 
#include <vector>
#include <iostream>
#include "ProtonEfficiency.h"
#include "TH2F.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

using namespace std;

ProtonEfficiency::ProtonEfficiency()
{
  //by default, a constant efficiency of 1.0
  constantEfficiencySet = true;
  constantEfficiency  = 1.0;  

  SetConstants();
}

ProtonEfficiency::ProtonEfficiency(TString filename)
{

  LoadTPCEff( filename );

  //by default, a constant efficiency of 1.0
  constantEfficiencySet = false;
  
  SetConstants();
}

ProtonEfficiency::ProtonEfficiency(TString filename1,TString filename2)
{

  LoadTPCEff( filename1 );
  LoadTOFEff( filename2 );

  //by default, a constant efficiency of 1.0
  constantEfficiencySet = false;
  
  SetConstants();

}


ProtonEfficiency::~ProtonEfficiency()
{
}


int ProtonEfficiency::SetConstants()
{
  _mass = 0.938272;
  _bin_edges = {-0.72 ,-0.64 ,-0.56 ,-0.48 ,-0.40 ,-0.32 ,-0.24 ,-0.16 ,-0.08 ,0.0,
		0.08  , 0.16 , 0.24, 0.32  , 0.40 ,0.48 , 0.56 , 0.64, 0.72 };
}

int ProtonEfficiency::LoadTOFEff(TString filename)
{
  cout << "LoadTOFEff() loading histogram..." << endl;
  TFile * f = new TFile(filename,"read");
  h_eff_tof = (TH2F*) f->Get("tpc_tofmatch");
  cout << "             histogram loaded." << endl;

  return 0;
}

int ProtonEfficiency::LoadTPCEff(TString filename)
{

  TFile * f = new TFile(filename,"read");

  cout << "LoadTPCEff() loading histograms and fits ..." << endl;
  for (int i=17;i<35;i++) 
    {
      _h1d_eff.push_back( (TH1D*) f->Get(TString::Format("projection_%i",i).Data()) );
      _tf1_eff.push_back( (TF1*) f->Get(TString::Format("f_%i",i).Data()) );
    }

  cout << "             hists and fits loaded" << endl;


  return 0;

}

double ProtonEfficiency::GetTPCEfficiency(double y, double pt)
{

  y = -1.0*y; // Switch Conventions

  if (y <= _bin_edges[0] || y >= _bin_edges[_bin_edges.size()-1] ) return -1; // if the rapidity is out of the rapidity range, return -1

  //Get the index for the correct histogram and fit
  int index =-1;
  for ( int iBin=0;iBin<_bin_edges.size()-1;iBin++)
    {
      if ( y > _bin_edges[iBin] && y <= _bin_edges[iBin+1] ) index=iBin;
    } 
  if (index<0) return -1;

  TH1D * h = _h1d_eff[index];
  TF1 * f = _tf1_eff[index];
  
  float xMin = h->GetBinLowEdge(h->FindFirstBinAbove(0.01));
  float xMax = 3.0;
  if ( pt <= xMin || pt >= xMax ) return -1;

  double eff = f->Eval( pt );

  if ( eff > 1.0 || eff < 0.01 ) return -1; //Set to 1% cut off
  
  return eff;
}

double ProtonEfficiency::GetTOFEfficiency(double y, double pt)
{
  return h_eff_tof->GetBinContent( h_eff_tof->FindBin(y,pt) );
}

double ProtonEfficiency::GetConstantEfficiency()
{
  return constantEfficiency;
}

void ProtonEfficiency::SetConstantEfficiency(double val)
{
  constantEfficiencySet = true;
  constantEfficiency     = val;
}

double ProtonEfficiency::GetEff(int charge, double pt,double pz)
{

  //an empty proton efficiency will return 1, no eff correction.  
  if (!constantEfficiencySet)
    {

      double y = Rapidity(pt,pz);
      double eff = GetTPCEfficiency(y,pt);

      if ( eff <= 0 || eff > 1.0 ) return -1;
      return eff;
    }

  return GetConstantEfficiency();
  
}

double ProtonEfficiency::GetEff(int charge, double pt,double pz,bool tofmatch)
{

  //an empty proton efficiency will return 1, no eff correction.  
  if (!constantEfficiencySet)
    {
      double y = Rapidity(pt,pz);       
      double eff = GetTPCEfficiency(y,pt);
      
      if (tofmatch)
	{
	  double tofEff = GetTOFEfficiency(y,pt);
	  eff = eff*tofEff;
	}
      
      if ( eff <= 0.01 || eff > 1.0 ) return -1;

      return eff;
    }

  return GetConstantEfficiency();
  
}


double ProtonEfficiency::Rapidity(double pt, double pz)
{
  Double_t energy = sqrt( pow(pt,2) + pow(pz,2) + pow(_mass,2) );
  Double_t rap = 0.5 * log( ( energy + pz) / (energy - pz) );
  rap = fabs(rap);
  
  return rap - 1.04;
}	  
