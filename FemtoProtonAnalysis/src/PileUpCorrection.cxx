#include <iostream>
#include <vector>
#include <utility>
#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "FlucContainer.h"
#include "PileUpCorrection.h"
#include "CumulantProfileContainer.h"

using std::cout;
using std::endl;
using std::vector;

PileUpCorrection::PileUpCorrection()
{
  int r = ORDER;

  for(int j=0;j<9;++j)
    {
      for(int i=0;i<r+1;++i)
	{
	  cbwcC[j][i]=0.;
	  C[j][i]=0.;
	  mu[j][i]=0.;
	  C_true[j][i]=0.;
	  cbwcC_true[j][i]=0.;
	  mu_true[j][i]=0.;
	  C_cor[j][i]=0.;
	  cbwcC_cor[j][i]=0.;
	  mu_cor[j][i]=0.;
	}
    }
}

PileUpCorrection::~PileUpCorrection()
{
}

CumulantProfileContainer * PileUpCorrection::LoadCumulant(TString filename,int iContainer)
{
  TFile * f = new TFile(filename,"read");
  
  //Setting C[][]
  CumulantProfileContainer * cpc = new CumulantProfileContainer(f,iContainer);  
  cpc->FactorialCumulantsToNormalCumulants();
  std::vector<TGraphErrors*> cumulantGraphs;
  cumulantGraphs.resize(5,NULL);
  TGraphErrors * wgr = cpc->GetWeightGraph();
  cumulantGraphs[1] = cpc->GetQ1Graph();
  cumulantGraphs[2] = cpc->GetQ2Graph();
  cumulantGraphs[3] = cpc->GetQ3Graph();
  cumulantGraphs[4] = cpc->GetQ4Graph();

  TH1F * mult = new TH1F("mult","mult",300,0,300);
  for (int i=0;i<wgr->GetN();i++)
    {
      double xMult = wgr->GetX()[i];
      double yVal = wgr->GetY()[i];
      mult->SetBinContent(mult->FindBin(xMult),yVal);
    }

  inputMult = mult;

  for ( int iMult=0;iMult<300;iMult++)
    {

      FlucContainer fluc;

      for (int iC=1;iC<5;iC++)
	{
	  C[iMult][iC] = cumulantGraphs[iC]->GetY()[iMult];
	  //	  C[iMult][iC] = 1.0;
	  fluc.SetCumulants(iC,C[iMult][iC]);
	}

      for (int iMu=1;iMu<5;iMu++)
	{
	  mu[iMult][iMu] = fluc.GetMomentsFromCumulants( iMu );
	}

    }
  return cpc;
}


CumulantProfileContainer * PileUpCorrection::LoadURQMDCumulant(TString filename)
{
  TFile * f = new TFile(filename,"read");
  
  inputMult = (TH1F*)f->Get("hMult_Measured");
  inputMult->SetName("inputMult");

  for ( int iMult=0;iMult<300;iMult++)
    {      
      for (int iMu=0;iMu<5;iMu++)
	{
	  mu[iMult][iMu] = 0;
	}

      if (iMult >= 199) continue;
      TH1D * hp_mult = (TH1D*)f->Get(Form("Measured%d",iMult));
      FlucContainer fluc;                                              
      LongDouble_t nevent = hp_mult->Integral();

      const int nxbin = hp_mult->GetXaxis()->GetNbins();
      for(int ixbin=0; ixbin<nxbin; ixbin++)
	{
	  const double nBinEvent = hp_mult->GetBinContent(ixbin+1);
	  if ( nBinEvent==0 ) continue;
	  const double m = hp_mult->GetBinCenter(ixbin+1);
	  for(int r=0; r<5;r++)
	    {
	      mu[iMult][r] += pow(m,r)*nBinEvent;
	    }
	}

      	for(int r=0; r<=ORDER; r++)
	  {
	    mu[iMult][r] /= nevent;
	    fluc.SetMoments(r,mu[iMult][r]);
	  }
	for (int iC=1;iC<5;iC++)
	  {
	    C[iMult][iC] = fluc.GetCumulant(iC);
	  }

	delete hp_mult;	
    }

  //  if ( hMult_true ) inputMult = hMult_true;

  CumulantProfileContainer * cpc = new CumulantProfileContainer(0,"UnCorr",(TH1F*)inputMult,C,300);
  return cpc;
}


CumulantProfileContainer * PileUpCorrection::GetURQMDTrueCumulant(TString filename)
{
  TFile * f = new TFile(filename,"read");
  
  TH1F * inputMult0 = (TH1F*)f->Get("hMult_True");
  inputMult->SetName("inputMult0");

  for ( int iMult=0;iMult<300;iMult++)
    {      
      for (int iMu=0;iMu<5;iMu++)
	{
	  mu_true[iMult][iMu] = 0;
	}

      if (iMult >= 199) continue;
      TH1D * hp_mult = (TH1D*)f->Get(Form("True%d",iMult));
      FlucContainer fluc;                                              
      LongDouble_t nevent = hp_mult->Integral();

      const int nxbin = hp_mult->GetXaxis()->GetNbins();
      for(int ixbin=0; ixbin<nxbin; ixbin++)
	{
	  const double nBinEvent = hp_mult->GetBinContent(ixbin+1);
	  if ( nBinEvent==0 ) continue;
	  const double m = hp_mult->GetBinCenter(ixbin+1);
	  for(int r=0; r<5;r++)
	    {
	      mu_true[iMult][r] += pow(m,r)*nBinEvent;
	    }
	}

      	for(int r=0; r<=ORDER; r++)
	  {
	    mu_true[iMult][r] /= nevent;
	    fluc.SetMoments(r,mu_true[iMult][r]);
	  }
	for (int iC=1;iC<5;iC++)
	  {
	    C_true[iMult][iC] = fluc.GetCumulant(iC);
	  }

	delete hp_mult;	
    }

  CumulantProfileContainer * cpc = new CumulantProfileContainer(0,"True",(TH1F*)inputMult0,C_true,300);
  return cpc;
}


bool PileUpCorrection::LoadPileUpHistograms(TString filename)
{
  TFile * f = new TFile(filename,"read");

  // Multiplicity dist, measured, true and reponse matrix
  /*
  hMult_true = (TH1F*)f->Get("FxtMult3_true"); // True
  hMult = (TH1F*)f->Get("FxtMult3");     // All
  hRM = (TH2F*)f->Get("FxtMult3_RM"); //// response matrix 
  */
  hMult_true = (TH1F*)f->Get("hMC_true"); // True
  hMult = (TH1F*)f->Get("hMC_all_final");     // All
  hRM = (TH2F*)f->Get("hRM_final"); //// response matrix 

  if ( hMult_true && hMult && hRM ) return true;
  else return false;
}

bool PileUpCorrection::LoadURQMDHistograms(TString filename)
{
  TFile * f = new TFile(filename,"read");

  // Multiplicity dist, measured, true and reponse matrix
  hMult_true = (TH1F*)f->Get("hMultTrue_scaled"); // True
  hMult = (TH1F*)f->Get("hMult_scaled");     // All
  hRM = (TH2F*)f->Get("hRM"); //// response matrix 

  if ( hMult_true && hMult && hRM ) return true;
  else return false;
}


void PileUpCorrection::CorrectionZero( LongDouble_t *mucor,LongDouble_t *musub, const LongDouble_t alpha, const LongDouble_t beta )
{
  for(int n=1; n<=4; n++){
    mucor[n] = 0;
  }
  mucor[1] = musub[1]/(1-alpha+2*beta);
  mucor[2] = (musub[2]-2*beta*mucor[1]*mucor[1])/(1-alpha+2*beta);
  mucor[3] = (musub[3]-6*beta*mucor[2]*mucor[1])/(1-alpha+2*beta);
  mucor[4] = (musub[4]-2*beta*(4*mucor[1]*mucor[3]+3*mucor[2]*mucor[2]))/(1-alpha+2*beta);
  
}

void PileUpCorrection::Correction( LongDouble_t *mucor,LongDouble_t *musub, LongDouble_t  *muzero, const LongDouble_t alpha, const LongDouble_t beta )
{
  for(int n=1; n<=4; n++){
    mucor[n] = 0;
  }
  mucor[1] = (musub[1]-beta*muzero[1])/(1-alpha+beta);
  mucor[2] = (musub[2]-beta*(muzero[2]+2*mucor[1]*muzero[1]))/(1-alpha+beta);
  mucor[3] = (musub[3]-beta*(muzero[3]+3*mucor[2]*muzero[1]+3*mucor[1]*muzero[2]))/(1-alpha+beta);
  mucor[4] = (musub[4]-beta*(muzero[4]+4*mucor[3]*muzero[1]+4*mucor[1]*muzero[3]+6*mucor[2]*muzero[2]))/(1-alpha+beta);
}

CumulantProfileContainer * PileUpCorrection::CorrectionForMultRange(int iName,int trgMult, int maxMult)
{

  const int r = ORDER;
  const int cutoff = 100;
  const int nbinx = hRM->GetXaxis()->GetNbins();
  const int nbiny = hRM->GetYaxis()->GetNbins();
	
  LongDouble_t muzero[r+1] = {};

  for (int iMult=0; iMult<maxMult; iMult++)
    {
      const LongDouble_t binevent = hMult_true->GetBinContent(iMult+1);
      //      if(binevent<cutoff) continue;
      
      //If below trig efficiency cut, set corrected cumulants to measured cumulants
      if (iMult<trgMult)
	{
	  for(int n=1; n<=r; n++) C_cor[iMult][n] = C[iMult][n];
	  continue;
	}
            
      LongDouble_t mu_pu[r+1] = {};
      //// Reconstruct moments for pile-up events      
      LongDouble_t weight_pu = 0;
      LongDouble_t weight_pu_true = 0;
      LongDouble_t weight_pu_all = 0;

      for(int ibin=0; ibin<nbinx; ibin++)
	{
	  for(int jbin=0; jbin<nbiny; jbin++)
	    {
	    if((ibin+jbin)!=iMult) continue;
	    weight_pu_all += hRM->GetBinContent(ibin+1,jbin+1)/hMult->GetBinContent(iMult+1);
	    // Skip if i==0 || j==0
	    if(ibin==iMult||jbin==iMult)
	      {
		weight_pu_true += hRM->GetBinContent(ibin+1,jbin+1)/hMult->GetBinContent(iMult+1);
		continue;
	      }
	    LongDouble_t mu_pu_sub[r+1] = {};
	    FlucContainer fpu;
	    // Superpose two collisions in terms of cumulants : mini-pileup
    	    for(int n=1; n<=r; n++) fpu.SetCumulants(n,C[ibin][n] + C[jbin][n]);
	    // Convert to moments to pile up those mini-pileups : pileup
	    for(int n=1; n<=r; n++)
	      {
		mu_pu_sub[n] = fpu.GetMomentsFromCumulants_v2(n);
		//		cout << " mu_pu_sub[" << n << "] = " << mu_pu_sub[n] << endl; 
	      }
	    // Weight for mini-pileup compared to whole events
	    const LongDouble_t weight = hRM->GetBinContent(ibin+1,jbin+1)/hMult->GetBinContent(iMult+1);
	    // Weight for pileup compared to whole events
	    weight_pu += weight;
	    // Propagate mini-pileup with proper weights 
	    for(int n=1; n<=r; n++) mu_pu[n] += weight*mu_pu_sub[n]; 
	  }
	}
      
      LongDouble_t musub[r+1] = {};
      for(int n=1; n<=r; n++) musub[n] = mu[iMult][n] - mu_pu[n];  // measured - pileup (i!=iMult||j!=iMult) 
      LongDouble_t mucor[r+1] = {};
      FlucContainer fcor;

      if(iMult==0)
	{
	  CorrectionZero(mucor,musub,weight_pu_all,weight_pu_true);
	  for(int n=1; n<=r; n++) muzero[n] = mucor[n];
	}

      if(iMult>0) Correction(mucor,musub,muzero,weight_pu_all,weight_pu_true);

      for(int n=1; n<=r; n++)
	{
	  mu_cor[iMult][n] = mucor[n]; 
	  fcor.SetMoments(n,mu_cor[iMult][n]);
	}

      //      cout <<"======================"<< endl;
      //      cout <<"Multiplicity #"<< iMult << endl;
      for(int n=1; n<=r; n++)
	{
	  C_cor[iMult][n] = fcor.GetCumulant(n);
	  //	  cout << n <<" "<< C[iMult][n] <<" "<< C_cor[iMult][n]  << endl;
	  C[iMult][n] = C_cor[iMult][n]; // replace measured cumulant to corrected cumulant
	}

    }

  inputMult = hMult_true;

  CumulantProfileContainer * cpc = new CumulantProfileContainer(iName,"Corr",(TH1F*)inputMult,C,maxMult);

  return cpc;
  
}


