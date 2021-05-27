//Class for proton efficiency
//Currently only has constant efficiency
//should be updated to have efficiency curves 
#include <vector>
#include <iostream>
#include <map>
#include "ProtonEfficiency.h"
#include "TH2F.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "InputParameterList.h"

using namespace std;

ProtonEfficiency::ProtonEfficiency()
{
  //by default, a constant efficiency of 1.0
  constantEfficiencySet = true;
  constantEfficiency  = 1.0;  

  SetConstants();
}

ProtonEfficiency::ProtonEfficiency(TString filename1,TString filename2,InputParameterList &pl)
{

  SetConstants();

  LoadTPCEff( filename1 );
  LoadTOFEff( filename2 );

  SetLabels( pl );

  //by default, a constant efficiency of 1.0
  constantEfficiencySet = false;
  
}

ProtonEfficiency::~ProtonEfficiency()
{
}


int ProtonEfficiency::SetConstants()
{
  _mass = 0.938272;

  for ( int i=0;i<93;i++)
    {
      _tpc_bin_edges.push_back( 0.02*i );
    }

  for ( int i=0;i<291;i++)
    {
      _tof_bin_edges.push_back( 0.005*i );
    }

  std::vector<TString> tofLabels = { "nhf_10", "nhf_12","nhf_15","m2_minus0.05","m2_plus0.05" };
  std::vector<TString> tpcLabels = { "nhf_10", "nhf_12","nhf_15" };

  tofSysLabels = tofLabels;
  tpcSysLabels = tpcLabels;

}

int ProtonEfficiency::LoadTOFEff(TString filename)
{

  TFile * f = new TFile(filename,"read");

  for ( auto &sysLabel : tofSysLabels )
    {
      
      for (int i=0;i<291;i++) 
	{
	  _tg_tof_eff[sysLabel].push_back( (TGraph*) NULL );

	  if ( i >= 15 ) _tg_tof_eff[sysLabel][i] = (TGraph*) f->Get(TString::Format("Interp/interpolate_Graph_%s_%i",sysLabel.Data(),i+1).Data());
	  //	  if ( !_tg_tof_eff[sysLabel][i] ) cout << "Cant find " << TString::Format("Interp/interpolate_Graph_%s_%i",sysLabel.Data(),i+1) << endl;
	}
    }

  return 0;
}

int ProtonEfficiency::LoadTPCEff(TString filename)
{

  TFile * f = new TFile(filename,"read");

  f->ls();

  for ( auto &sysLabel : tpcSysLabels )
    {

      for (int i=0;i<92;i++) 
	{
	  _tg_tpc_eff[sysLabel].push_back( (TGraphAsymmErrors*) NULL );
	  _tf1_tpc_eff[sysLabel].push_back((TF1*) NULL );
	  _xMin[sysLabel].push_back( 0 );
	  if ( i >= 2 )_tg_tpc_eff[sysLabel][i] = (TGraphAsymmErrors*) f->Get(TString::Format("Ratio_Sys_%s_slice_bin%i",sysLabel.Data(),i+1).Data());
	  _tf1_tpc_eff[sysLabel].push_back( (TF1*) NULL );
	  if ( i >= 2 ) _tf1_tpc_eff[sysLabel][i] = (TF1*) f->Get(TString::Format("f_Sys_%s_slice_bin%i",sysLabel.Data(),i+1).Data());
	  
	  if ( !_tf1_tpc_eff[sysLabel][i] )
	    {
	      //	      cout << "Cant find " << TString::Format("f_Sys_%s_slice_bin%i",sysLabel.Data(),i+1) << endl;
	      continue;
	    }
	  //Get first point above 1%
	  for (int ipt=0;ipt<_tg_tpc_eff[sysLabel][i]->GetN();ipt++)
	    {
	      if ( _tg_tpc_eff[sysLabel][i]->GetX()[ipt] > 0.01 ) {
		_xMin[sysLabel][i] = _tg_tpc_eff[sysLabel][i]->GetX()[ipt];
		break;
	      }
	    }
	}
    }

  return 0;

}

double ProtonEfficiency::GetTPCEfficiency(double y, double pt)
{
    
  if (y <= _tpc_bin_edges[2] || y >= _tpc_bin_edges[_tpc_bin_edges.size()-1] ) return -1; // if the rapidity is out of the rapidity range, return -1

  //Get the index for the correct histogram and fit
  int index =-1;
  for ( int iBin=2;iBin<_tpc_bin_edges.size()-1;iBin++)
    {
      if ( y > _tpc_bin_edges[iBin] && y <= _tpc_bin_edges[iBin+1] ) index=iBin;
    } 

  if (index<0) return -1;
  TGraphAsymmErrors * gr = _tg_tpc_eff[tpcLabel][index];

  TF1 * f = _tf1_tpc_eff[tpcLabel][index];
 
  float xMin = _xMin[tpcLabel][index];
  float xMax = 2.2;

  if ( pt <= xMin || pt >= xMax ) return -1;
  double eff = f->Eval( pt );
  if ( eff > 1.0 || eff < 0.01 ) return -1; //Set to 1% cut off
  
  return eff;
}

double ProtonEfficiency::GetTOFEfficiency(double eta, double pt)
{

  if (eta <= _tof_bin_edges[15] || eta >= _tof_bin_edges[_tof_bin_edges.size()-1] ) return -1; // if the rapidity is out of the rapidity range, return -1

  //Get the index for the correct histogram and fit
  int index =-1;
  for ( int iBin=15;iBin<_tof_bin_edges.size()-1;iBin++)
    {
      if ( eta > _tof_bin_edges[iBin] && eta <= _tof_bin_edges[iBin+1] ) index=iBin;
    } 
  if (index<0) return -1;

  TGraph * gr = _tg_tof_eff[tofLabel][index];

  float xMin = 0.5;
  float xMax = 2.2;
  if ( pt <= xMin || pt >= xMax ) return -1;
  double eff = gr->Eval( pt );
  if ( eff > 1.0 || eff < 0.01 ) return -1; //Set to 1% cut off
  
  return eff;
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

double ProtonEfficiency::GetEff(int charge, double pt,double pz,bool tofmatch)
{

  //an empty proton efficiency will return 1, no eff correction.  
  if (!constantEfficiencySet)
    {
      double eff=0;
      double y = Rapidity(pt,pz);       
      double tpcEff = GetTPCEfficiency(y,pt);
      
      if (tofmatch)
	{
	  double eta = Pseudorapidity(pt,pz);       
	  double tofEff = GetTOFEfficiency(eta,pt);
	  eff = tpcEff*tofEff;
	}
      else eff = tpcEff;

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
  return rap;
}	  

double ProtonEfficiency::Pseudorapidity(double pt, double pz)
{
  Double_t p = sqrt(pt*pt + pz*pz);
  Double_t eta = 0.5*TMath::Log( ( p + pz )/( p - pz ) );
  eta = fabs(eta);
  return eta;
}	  


void ProtonEfficiency::SetLabels(InputParameterList &pl)
{
  if ( pl.NHitsFitMin() == 10 ) 
    {
      tpcLabel = tpcSysLabels[0];

      if ( pl.Mass2Low() <= 0.55 ) tofLabel = tofSysLabels[3];
      if ( pl.Mass2Low() >= 0.65 ) tofLabel = tofSysLabels[4];
      else tofLabel = tofSysLabels[0];
    }  
  if ( pl.NHitsFitMin() == 12 )
    {
      tofLabel = tofSysLabels[1];
      tpcLabel = tofSysLabels[2];
    }

  if ( pl.NHitsFitMin() == 15 )
    {
      tofLabel = tofSysLabels[2];
      tpcLabel = tpcSysLabels[2];
    }
}

void ProtonEfficiency::UnitTest(TH2D * tpcHist, TH2D * tofHist)
{

  cout << "******************************************" << endl;
  cout << "******************************************" << endl;
  cout << "**************UNIT TESTING****************" << endl;
  cout << "******************************************" << endl;
  cout << "******************************************" << endl;

  TH2D * tpc_check = (TH2D*) tpcHist->Clone();
  tpc_check->Reset();
  tpc_check->SetName("TPCratio");

  cout << "LOAD TPC FILE " << endl;

  cout << "TPC LABEL = " << tpcLabel << endl;
  cout << "TOF LABEL = " << tofLabel << endl;
  
  for ( int iX=0;iX<tpc_check->GetNbinsX();iX++)
    {
      for ( int iY=0;iY<tpc_check->GetNbinsY();iY++)
  	{
  
	  double x = tpc_check->GetXaxis()->GetBinCenter(iX);
	  double y = tpc_check->GetYaxis()->GetBinCenter(iY);
	  double binContent = tpcHist->GetBinContent(iX,iY);
	  
	  double eff = GetTPCEfficiency( x ,y );

	  if ( eff > 0 ) tpc_check->SetBinContent(iX,iY, binContent/eff);
	}
    }

  cout << "LOAD TOF FILE " << endl;
  
  TH2D * tof_check = (TH2D*) tofHist->Clone();
  tof_check->Reset();
  tof_check->SetName("TOFratio");

  for ( int iX=0;iX<tof_check->GetNbinsX();iX++)
    {
      for ( int iY=0;iY<tof_check->GetNbinsY();iY++)
	{
	  double x = tof_check->GetXaxis()->GetBinCenter(iX);
	  double y = tof_check->GetYaxis()->GetBinCenter(iY);
	  double x_low = tof_check->GetXaxis()->GetBinLowEdge(iX);
	  double x_high = tof_check->GetXaxis()->GetBinLowEdge(iX) + tof_check->GetXaxis()->GetBinWidth(iX);
	  double binContent = tofHist->GetBinContent(iX,iY);

	  cout << "eta = [" << x_low << "," << x_high << "]" << endl; 
	  double eff = GetTOFEfficiency( x ,y );
	  
	  cout << endl;
	  if ( eff > 0 ) tof_check->SetBinContent(iX,iY, binContent/eff);
	}
    }

  TFile * UnitFile = new TFile("Eff_UnitTest.root","recreate");
  UnitFile->cd();
  tof_check->Write();
  tpc_check->Write();
  UnitFile->Close();

}
