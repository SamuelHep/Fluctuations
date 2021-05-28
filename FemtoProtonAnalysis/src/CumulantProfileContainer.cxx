
/*
CumulantProfileContainer is the main analysis class.
What the class does:
1. Takes bivariate moments (m_r_s) from an event and calculates
   all the needed variables i.e: m_1_1, m_1_1*m_1_1 ... m_1_1*m_2_1*m_2_1.
2. bivariate moment variables are put in profiles
3. After running through every event, profiles are converted from moments to cumulants
4. Cumulants are rebinned with CBWC method and cumulant ratios are calcuted

The class is designed to have steps 1 and 2 run through condor, while 3 and 4 are run locally.

The code to calculate m_r_s is in src/MomentFunctions.cxx
 */

#include <iostream>
#include <vector>
#include <utility>
#include <map>

#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "CumulantFunctions.h"

#include "CumulantProfileContainer.h"

ClassImp(CumulantProfileContainer)

using namespace std;

//Basic Constructor
CumulantProfileContainer::CumulantProfileContainer(int iName)
{
  
  InitProfiles(iName);
  InitGraphs(iName); // initialized _gr[""] 

}

//Takes an array of cumulants and fills the CumulantProfileContainer Graphs (C1, C2, C3, ... )
//This is used after in the pileup correction code
CumulantProfileContainer::CumulantProfileContainer(int iName,TString label,TH1F * multHist,LongDouble_t Cumulants[300][7],int MaxMult)
{

  InitGraphs(iName); // initialized _gr[""] 

  for (int iMult=0;iMult<MaxMult;iMult++)
    {

      long double C1 = Cumulants[iMult][1];
      long double C2 = Cumulants[iMult][2];
      long double C3 = Cumulants[iMult][3];
      long double C4 = Cumulants[iMult][4];
      long double C5 = Cumulants[iMult][5];
      long double C6 = Cumulants[iMult][6];

      long double K1 = C1;
      long double K2 = C2 - K1;
      long double K3 = C3 - 3*C2 +2*C1;
      long double K4 = -6*C1 + 11*C2 -6*C3 + C4;
      long double K5 = -10*C4 + 35*C3 -50*C2 + 24*C1 + C5;
      long double K6 = -120*C1 + 274*C2 - 225*C3 + 85*C4 - 15*C5 + C6;

      _gr["N"]->SetPoint(_gr["N"]->GetN(),iMult,multHist->GetBinContent(multHist->FindBin(iMult)));
      _gr["C1"]->SetPoint(_gr["C1"]->GetN(),iMult,C1);
      _gr["C2"]->SetPoint(_gr["C2"]->GetN(),iMult,C2);
      _gr["C3"]->SetPoint(_gr["C3"]->GetN(),iMult,C3);
      _gr["C4"]->SetPoint(_gr["C4"]->GetN(),iMult,C4);
      _gr["C5"]->SetPoint(_gr["C5"]->GetN(),iMult,C5);
      _gr["C6"]->SetPoint(_gr["C6"]->GetN(),iMult,C6);

      _gr["K1"]->SetPoint(_gr["K1"]->GetN(),iMult,K1);
      _gr["K2"]->SetPoint(_gr["K2"]->GetN(),iMult,K2);
      _gr["K3"]->SetPoint(_gr["K3"]->GetN(),iMult,K3);
      _gr["K4"]->SetPoint(_gr["K4"]->GetN(),iMult,K4);
      _gr["K5"]->SetPoint(_gr["K5"]->GetN(),iMult,K5);
      _gr["K6"]->SetPoint(_gr["K6"]->GetN(),iMult,K6);

    }
  
}

//Inits the class from a TFile with filed profiles
CumulantProfileContainer::CumulantProfileContainer(TFile *f,int iName)
{
  InitProfiles(iName,f); //Load profiles
  InitGraphs(iName); // initialized _gr[""] 
}

//Destructor
CumulantProfileContainer::~CumulantProfileContainer()
{
}

//If you want to change the names of your graphs C1 -> C1_uncor
void CumulantProfileContainer::AmendGraphSuffix( TString suffix ) 
{
  for ( auto &gr_map : _gr ) 
    {
      gr_map.second->SetName( TString::Format("%s%s",gr_map.second->GetName(), suffix.Data() ) );
    }
}

//Initializes all the graphs and sets their names
void CumulantProfileContainer::InitGraphs(int iName)
{

  TString ext = ( iName == kPrimary ) ? "" : TString::Format("_index%i",iName); // Add extension if bootsrap

  _gr["N"] = new TGraphErrors(); 
  _gr["C1"] = new TGraphErrors();
  _gr["C2"] = new TGraphErrors();
  _gr["C3"] = new TGraphErrors();
  _gr["C4"] = new TGraphErrors();
  _gr["C5"] = new TGraphErrors();
  _gr["C6"] = new TGraphErrors();

  _gr["C2_C1"] = new TGraphErrors();
  _gr["C3_C2"] = new TGraphErrors();
  _gr["C4_C2"] = new TGraphErrors();
  _gr["C5_C1"] = new TGraphErrors();
  _gr["C6_C2"] = new TGraphErrors();

  _gr["C1_cbwc"] = new TGraphErrors();
  _gr["C2_cbwc"] = new TGraphErrors();
  _gr["C3_cbwc"] = new TGraphErrors();
  _gr["C4_cbwc"] = new TGraphErrors();
  _gr["C5_cbwc"] = new TGraphErrors();
  _gr["C6_cbwc"] = new TGraphErrors();

  _gr["C2_C1"] = new TGraphErrors();
  _gr["C3_C2"] = new TGraphErrors();
  _gr["C4_C2"] = new TGraphErrors();
  _gr["C5_C1"] = new TGraphErrors();
  _gr["C6_C2"] = new TGraphErrors();

  _gr["C2_C1_cbwc"] = new TGraphErrors();
  _gr["C3_C2_cbwc"] = new TGraphErrors();
  _gr["C4_C2_cbwc"] = new TGraphErrors();
  _gr["C5_C1_cbwc"] = new TGraphErrors();
  _gr["C6_C2_cbwc"] = new TGraphErrors();

  _gr["K1"] = new TGraphErrors();
  _gr["K2"] = new TGraphErrors();
  _gr["K3"] = new TGraphErrors();
  _gr["K4"] = new TGraphErrors();
  _gr["K5"] = new TGraphErrors();
  _gr["K6"] = new TGraphErrors();

  _gr["K2_K1"] = new TGraphErrors();
  _gr["K3_K1"] = new TGraphErrors();
  _gr["K4_K1"] = new TGraphErrors();
  _gr["K5_K1"] = new TGraphErrors();
  _gr["K6_K1"] = new TGraphErrors();

  _gr["K1_cbwc"] = new TGraphErrors();
  _gr["K2_cbwc"] = new TGraphErrors();
  _gr["K3_cbwc"] = new TGraphErrors();
  _gr["K4_cbwc"] = new TGraphErrors();
  _gr["K5_cbwc"] = new TGraphErrors();
  _gr["K6_cbwc"] = new TGraphErrors();

  _gr["K2_K1_cbwc"] = new TGraphErrors();
  _gr["K3_K1_cbwc"] = new TGraphErrors();
  _gr["K4_K1_cbwc"] = new TGraphErrors();
  _gr["K5_K1_cbwc"] = new TGraphErrors();
  _gr["K6_K1_cbwc"] = new TGraphErrors();

  //Set All Titles
  for ( auto &gr_m : _gr )
    {
      gr_m.second->SetName( TString::Format("%s%s",gr_m.first.Data(),ext.Data()) );
    }

}

//Set all the profiles to NULL and then load or make new profiles
void CumulantProfileContainer::InitProfiles(int iName, TFile * file)
{

  _profile["m11"] = NULL;
  _profile["m11_2"] = NULL;
  _profile["m21_2"] = NULL;
  _profile["m21"] = NULL;
  _profile["m22"] = NULL;
  _profile["m11_3"] = NULL;
  _profile["m11_m21"] = NULL;
  _profile["m11_m22"] = NULL;
  _profile["m31"] = NULL;
  _profile["m32"] = NULL;
  _profile["m33"] = NULL;
  _profile["m11_4"] = NULL;
  _profile["m11_2_m21"] = NULL;
  _profile["m11_2_m22"] = NULL;
  _profile["m11_m31"] = NULL;
  _profile["m21_2"] = NULL;
  _profile["m22_2"] = NULL;
  _profile["m11_m32"] = NULL;
  _profile["m11_m33"] = NULL;
  _profile["m21_m22"] = NULL;
  _profile["m41"] = NULL;
  _profile["m42"] = NULL;
  _profile["m43"] = NULL;
  _profile["m44"] = NULL;
  
  _profile["m51"] = NULL;
  _profile["m52"] = NULL;
  _profile["m53"] = NULL;
  _profile["m54"] = NULL;
  _profile["m55"] = NULL;
  _profile["m11_m41"] = NULL;
  _profile["m11_m42"] = NULL;
  _profile["m11_m43"] = NULL;
  _profile["m11_m44"] = NULL;
  _profile["m21_m31"] = NULL;
  _profile["m21_m32"] = NULL;
  _profile["m21_m33"] = NULL;
  _profile["m22_m31"] = NULL;
  _profile["m22_m32"] = NULL;
  _profile["m22_m33"] = NULL;
  
  _profile["m11_2_m31"] = NULL;
  _profile["m11_2_m32"] = NULL;
  _profile["m11_2_m33"] = NULL;
  _profile["m22_2_m11"] = NULL;
  _profile["m21_2_m11"] = NULL;
  _profile["m11_m21_m22"] = NULL;
  
  _profile["m11_3_m21"] = NULL;
  _profile["m11_3_m22"] = NULL;
  _profile["m11_5"] = NULL;
  
  _profile["m11_6"] = NULL;
  
  _profile["m11_4_m21"] = NULL;
  _profile["m11_4_m22"] = NULL;
  
  _profile["m11_3_m31"] = NULL;
  _profile["m11_3_m32"] = NULL;
  _profile["m11_3_m33"] = NULL;
  
  _profile["m11_2_m22_m21"] = NULL;
  _profile["m11_2_m21_2"] = NULL;
  _profile["m11_2_m22_2"] = NULL;
  
  _profile["m21_3"] = NULL;
  _profile["m22_3"] = NULL;
  _profile["m11_2_m41"] = NULL;
  
  _profile["m11_2_m42"] = NULL;
  _profile["m11_2_m43"] = NULL;
  _profile["m11_2_m44"] = NULL;
  
  _profile["m21_2_m22"] = NULL;
  _profile["m22_2_m21"] = NULL;
  
  _profile["m11_m21_m31"] = NULL;
  _profile["m11_m21_m32"] = NULL;
  _profile["m11_m21_m33"] = NULL;
  
  _profile["m11_m22_m31"] = NULL;
  _profile["m11_m22_m32"] = NULL;
  _profile["m11_m22_m33"] = NULL;
  
  _profile["m11_m51"] = NULL;
  _profile["m11_m52"] = NULL;
  _profile["m11_m53"] = NULL;
  _profile["m11_m54"] = NULL;
  _profile["m11_m55"] = NULL;
  
  _profile["m21_m41"] = NULL;
  _profile["m21_m42"] = NULL;
  _profile["m21_m43"] = NULL;
  _profile["m21_m44"] = NULL;
  
  _profile["m22_m41"] = NULL;
  _profile["m22_m42"] = NULL;
  _profile["m22_m43"] = NULL;
  _profile["m22_m44"] = NULL;
  
  _profile["m31_2"] = NULL;
  _profile["m31_m32"] = NULL;
  _profile["m31_m33"] = NULL;
  _profile["m32_2"] = NULL;
  _profile["m32_m33"] = NULL;
  _profile["m33_2"] = NULL;
  
  _profile["m61"] = NULL;
  _profile["m62"] = NULL;
  _profile["m63"] = NULL;
  _profile["m64"] = NULL;
  _profile["m65"] = NULL;
  _profile["m66"] = NULL;

  TString ext = ( iName == kPrimary ) ? "" : TString::Format("_INDEX%i",iName);
  
  if ( file == NULL )   //MAKE NEW PROFILES
    {
      for ( auto &m : _profile )
	{
	  m.second = new TProfile( TString::Format("Profile_%s%s", m.first.Data(),ext.Data()) , "", 300, -0.5, 299.5);
	}
    }
  else 
    {
      for ( auto &m : _profile ) // LOAD PROFILES
	{
	  m.second = (TProfile*) file->Get(TString::Format("Profile_%s%s", m.first.Data(),ext.Data()));
	}
    }
  
}

//Calculate the m_r_s variables to be averaged 
void CumulantProfileContainer::FillBiVariateMoments(std::vector<std::vector<long double>> & m_r_s)
{

  m["m11"] = m_r_s[1][1]; // m(1,1)
  m["m11_2"] = m_r_s[1][1]*m_r_s[1][1]; //  m(1,1)^2 
  m["m21_2"] = m_r_s[2][1]*m_r_s[2][1]; //  m(2,1)^2 
  m["m21"] = m_r_s[2][1]; //  m(2,1)^2 
  m["m22"] = m_r_s[2][2]; 
  m["m11_3"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1];
  m["m11_m21"] = m_r_s[1][1]*m_r_s[2][1];
  m["m11_m22"] = m_r_s[1][1]*m_r_s[2][2];
  m["m31"] = m_r_s[3][1];
  m["m32"] = m_r_s[3][2];
  m["m33"] = m_r_s[3][3];
  m["m11_4"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1];
  m["m11_2_m21"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][1];
  m["m11_2_m22"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][2];
  m["m11_m31"] = m_r_s[1][1]*m_r_s[3][1];
  m["m21_2"] = m_r_s[2][1]*m_r_s[2][1];
  m["m22_2"] = m_r_s[2][2]*m_r_s[2][2];
  m["m11_m32"] = m_r_s[1][1]*m_r_s[3][2];
  m["m11_m33"] = m_r_s[1][1]*m_r_s[3][3];
  m["m21_m22"] = m_r_s[2][1]*m_r_s[2][2];
  m["m41"] = m_r_s[4][1];
  m["m42"] = m_r_s[4][2];
  m["m43"] = m_r_s[4][3];
  m["m44"] = m_r_s[4][4];

  m["m51"] = m_r_s[5][1];
  m["m52"] = m_r_s[5][2];
  m["m53"] = m_r_s[5][3];
  m["m54"] = m_r_s[5][4];
  m["m55"] = m_r_s[5][5];
  m["m11_m41"] = m_r_s[1][1]*m_r_s[4][1];
  m["m11_m42"] = m_r_s[1][1]*m_r_s[4][2];
  m["m11_m43"] = m_r_s[1][1]*m_r_s[4][3];
  m["m11_m44"] = m_r_s[1][1]*m_r_s[4][4];
  m["m21_m31"] = m_r_s[2][1]*m_r_s[3][1];
  m["m21_m32"] = m_r_s[2][1]*m_r_s[3][2];
  m["m21_m33"] = m_r_s[2][1]*m_r_s[3][3];
  m["m22_m31"] = m_r_s[2][2]*m_r_s[3][1];
  m["m22_m32"] = m_r_s[2][2]*m_r_s[3][2];
  m["m22_m33"] = m_r_s[2][2]*m_r_s[3][3];

  m["m11_2_m31"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[3][1];
  m["m11_2_m32"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[3][2];
  m["m11_2_m33"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[3][3];
  m["m22_2_m11"] = m_r_s[2][2]*m_r_s[2][2]*m_r_s[1][1];
  m["m21_2_m11"] = m_r_s[2][1]*m_r_s[2][1]*m_r_s[1][1];
  m["m11_m21_m22"] = m_r_s[1][1]*m_r_s[2][1]*m_r_s[2][2];

  m["m11_3_m21"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][1];
  m["m11_3_m22"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][2];
  m["m11_5"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1];

  m["m11_6"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1];
  
  m["m11_4_m21"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][1];
  m["m11_4_m22"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][2];

  m["m11_3_m31"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[3][1];
  m["m11_3_m32"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[3][2];
  m["m11_3_m33"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[3][3];

  m["m11_2_m22_m21"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][2]*m_r_s[2][1];
  m["m11_2_m21_2"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][1]*m_r_s[2][1];
  m["m11_2_m22_2"] =m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][2]*m_r_s[2][2];
  
  m["m21_3"] = m_r_s[2][1]*m_r_s[2][1]*m_r_s[2][1];
  m["m22_3"] = m_r_s[2][2]*m_r_s[2][2]*m_r_s[2][2];
  m["m11_2_m41"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[4][1];

  m["m11_2_m42"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[4][2];
  m["m11_2_m43"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[4][3];
  m["m11_2_m44"] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[4][4];

  m["m21_2_m22"] = m_r_s[2][1]*m_r_s[2][1]*m_r_s[2][2];
  m["m22_2_m21"] = m_r_s[2][2]*m_r_s[2][2]*m_r_s[2][1];

  m["m11_m21_m31"] = m_r_s[1][1]*m_r_s[2][1]*m_r_s[3][1]; 
  m["m11_m21_m32"] = m_r_s[1][1]*m_r_s[2][1]*m_r_s[3][2];
  m["m11_m21_m33"] = m_r_s[1][1]*m_r_s[2][1]*m_r_s[3][3];

  m["m11_m22_m31"] = m_r_s[1][1]*m_r_s[2][2]*m_r_s[3][1];
  m["m11_m22_m32"] = m_r_s[1][1]*m_r_s[2][2]*m_r_s[3][2];
  m["m11_m22_m33"] = m_r_s[1][1]*m_r_s[2][2]*m_r_s[3][3];

  m["m11_m51"] = m_r_s[1][1]*m_r_s[5][1];
  m["m11_m52"] = m_r_s[1][1]*m_r_s[5][2];
  m["m11_m53"] = m_r_s[1][1]*m_r_s[5][3];
  m["m11_m54"] = m_r_s[1][1]*m_r_s[5][4];
  m["m11_m55"] = m_r_s[1][1]*m_r_s[5][5];

  m["m21_m41"] = m_r_s[2][1]*m_r_s[4][1];
  m["m21_m42"] = m_r_s[2][1]*m_r_s[4][2];
  m["m21_m43"] = m_r_s[2][1]*m_r_s[4][3];
  m["m21_m44"] = m_r_s[2][1]*m_r_s[4][4];

  m["m22_m41"] = m_r_s[2][2]*m_r_s[4][1];
  m["m22_m42"] = m_r_s[2][2]*m_r_s[4][2];
  m["m22_m43"] = m_r_s[2][2]*m_r_s[4][3];
  m["m22_m44"] = m_r_s[2][2]*m_r_s[4][4];
  
  m["m31_2"] = m_r_s[3][1]*m_r_s[3][1];
  m["m31_m32"] = m_r_s[3][1]*m_r_s[3][2];
  m["m31_m33"] = m_r_s[3][1]*m_r_s[3][3];
  m["m32_2"] = m_r_s[3][2]*m_r_s[3][2];
  m["m32_m33"] = m_r_s[3][2]*m_r_s[3][3];
  m["m33_2"] = m_r_s[3][3]*m_r_s[3][3];
  
  m["m61"] = m_r_s[6][1];
  m["m62"] = m_r_s[6][2];
  m["m63"] = m_r_s[6][3];
  m["m64"] = m_r_s[6][4];
  m["m65"] = m_r_s[6][5];
  m["m66"] = m_r_s[6][6];

}

//Fill m_r_s for a give centrality bin
void CumulantProfileContainer::FillProfile(int cent, std::vector<std::vector<long double>> m_r_s)
{

  FillBiVariateMoments(m_r_s);
  for ( auto &prof_map : _profile )
    {
      prof_map.second->Fill( cent, m[ prof_map.first ] );
    }

  m.clear();

}

//Bivariate moments to bivariate cumulants 
void CumulantProfileContainer::InitCumulantMap(int cent)
{
  
  q.clear();
  for ( auto & val : m ) { val.second = -999; }

  for ( auto &_key_profile : _profile )
    {
      TString key = _key_profile.first;
      TProfile * profile = _key_profile.second;
      long double val = profile->GetBinContent( profile->FindBin(cent) );
      m[ key ] = val;
    }

  q["q11"]         = m["m11"];
  q["q11_2"]       = q_ab(m["m11_2"], m["m11"], m["m11"]);
  q["q21"]         = m["m21"];
  q["q22"]         = m["m22"]; 
  q["q11_3"]       = q_abc( m["m11_3"], m["m11_2"], m["m11_2"], m["m11_2"], m["m11"], m["m11"], m["m11"] );
  q["q11_q21"]     = q_ab( m["m11_m21"], m["m11"], m["m21"]);
  q["q11_q22"]     = q_ab( m["m11_m22"], m["m11"], m["m22"]);
  q["q31"]         = m["m31"];
  q["q32"]         = m["m32"];
  q["q33"]         = m["m33"];
  q["q11_4"]       = q_abcd( m["m11_4"],
			     m["m11_3"], m["m11_3"], m["m11_3"], m["m11_3"],
			     m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"],
			     m["m11"], m["m11"], m["m11"], m["m11"] );

  long double q11_4 = m["m11_4"] - 4*m["m11_3"]*m["m11"] -3*m["m11_2"]*m["m11_2"] + //q(1,1)^4
    12*m["m11_2"]*m["m11"]*m["m11"] -6*m["m11"]*m["m11"]*m["m11"]*m["m11"];
  

  q["q11_2_q21"]   = q_abc( m["m11_2_m21"], m["m11_2"], m["m11_m21"], m["m11_m21"], m["m11"], m["m11"], m["m21"] );
  q["q11_2_q22"]   = q_abc( m["m11_2_m22"], m["m11_2"], m["m11_m22"], m["m11_m22"], m["m11"], m["m11"], m["m22"] );
  q["q11_q31"]     = q_ab( m["m11_m31"], m["m11"], m["m31"]);
  q["q21_2"]       = q_ab( m["m21_2"], m["m21"], m["m21"]);
  q["q22_2"]       = q_ab( m["m22_2"], m["m22"], m["m22"]);
  q["q11_q32"]     = q_ab( m["m11_m32"], m["m11"], m["m32"]);
  q["q11_q33"]     = q_ab( m["m11_m33"], m["m11"], m["m33"]);
  q["q21_q22"]     = q_ab( m["m21_m22"], m["m21"], m["m22"]);
  q["q41"]         = m["m41"];
  q["q42"]         = m["m42"];
  q["q43"]         = m["m43"];
  q["q44"]         = m["m44"]; // q(4,4)

  q["q51"]         = m["m51"];
  q["q52"]         = m["m52"];
  q["q53"]         = m["m53"];
  q["q54"]         = m["m54"];
  q["q55"]         = m["m55"];

  q["q11_q41"]     = q_ab( m["m11_m41"], m["m11"], m["m41"] );
  q["q11_q42"]     = q_ab( m["m11_m42"], m["m11"], m["m42"] );
  q["q11_q43"]     = q_ab( m["m11_m43"], m["m11"], m["m43"] );
  q["q11_q44"]     = q_ab( m["m11_m44"], m["m11"], m["m44"] );
  q["q21_q31"]     = q_ab( m["m21_m31"], m["m21"], m["m31"] );
  q["q21_q32"]     = q_ab( m["m21_m32"], m["m21"], m["m32"] );
  q["q21_q33"]     = q_ab( m["m21_m33"], m["m21"], m["m33"] );
  q["q22_q31"]     = q_ab( m["m22_m31"], m["m22"], m["m31"] );
  q["q22_q32"]     = q_ab( m["m22_m32"], m["m22"], m["m32"] );
  q["q22_q33"]     = q_ab( m["m22_m33"], m["m22"], m["m33"] );

  q["q11_2_q31"]   = q_abc( m["m11_2_m31"], m["m11_2"], m["m11_m31"], m["m11_m31"], m["m11"], m["m11"], m["m31"] );
  q["q11_2_q32"]   = q_abc( m["m11_2_m32"], m["m11_2"], m["m11_m32"], m["m11_m32"], m["m11"], m["m11"], m["m32"] );
  q["q11_2_q33"]   = q_abc( m["m11_2_m33"], m["m11_2"], m["m11_m33"], m["m11_m33"], m["m11"], m["m11"], m["m33"] );
  q["q22_2_q11"]   = q_abc( m["m22_2_m11"], m["m22_2"], m["m11_m22"], m["m11_m22"], m["m22"], m["m22"], m["m11"] );
  q["q21_2_q11"]   = q_abc( m["m21_2_m11"], m["m21_2"], m["m11_m21"], m["m11_m21"], m["m21"], m["m21"], m["m11"] );
  q["q11_q21_q22"] = q_abc( m["m11_m21_m22"], m["m11_m21"], m["m11_m22"], m["m21_m22"], m["m11"], m["m21"], m["m22"] );

  q["q11_3_q21"]   = q_abcd( m["m11_3_m21"],
			     m["m11_3"], m["m11_2_m21"], m["m11_2_m21"], m["m11_2_m21"],
			     m["m11_2"], m["m11_2"], m["m11_m21"], m["m11_2"], m["m11_m21"], m["m11_m21"],
			     m["m11"], m["m11"], m["m11"], m["m21"] );
  q["q11_3_q22"]   = q_abcd( m["m11_3_m22"],
			     m["m11_3"], m["m11_2_m22"], m["m11_2_m22"], m["m11_2_m22"],
			     m["m11_2"], m["m11_2"], m["m11_m22"], m["m11_2"], m["m11_m22"], m["m11_m22"],
			     m["m11"], m["m11"], m["m11"], m["m22"] );


  q["q11_5"]       = q_abcde( m["m11_5"],
			      m["m11_4"], m["m11_4"], m["m11_4"], m["m11_4"], m["m11_4"],
			      m["m11_3"], m["m11_3"], m["m11_3"], m["m11_3"], m["m11_3"],
			      m["m11_3"], m["m11_3"], m["m11_3"], m["m11_3"], m["m11_3"],
			      m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"],
			      m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"],
			      m["m11"], m["m11"], m["m11"], m["m11"], m["m11"] );
  
  //Values for 6th order cumulant

  q["q11_4_q21"] = q_abcde( m["m11_4_m21"],
			    m["m11_4"], m["m11_3_m21"], m["m11_3_m21"], m["m11_3_m21"], m["m11_3_m21"],
			    m["m11_3"], m["m11_3"], m["m11_3"], m["m11_3"], m["m11_2_m21"],
			    m["m11_2_m21"], m["m11_2_m21"], m["m11_2_m21"], m["m11_2_m21"], m["m11_2_m21"],
			    m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"],
			    m["m11_2"], m["m11_m21"], m["m11_m21"], m["m11_m21"], m["m11_m21"],
			    m["m11"], m["m11"], m["m11"], m["m11"], m["m21"] );

  q["q11_4_q22"] = q_abcde( m["m11_4_m22"],
			    m["m11_4"], m["m11_3_m22"], m["m11_3_m22"], m["m11_3_m22"], m["m11_3_m22"],
			    m["m11_3"], m["m11_3"], m["m11_3"], m["m11_3"], m["m11_2_m22"],
			    m["m11_2_m22"], m["m11_2_m22"], m["m11_2_m22"], m["m11_2_m22"], m["m11_2_m22"],
			    m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"], m["m11_2"],
			    m["m11_2"], m["m11_m22"], m["m11_m22"], m["m11_m22"], m["m11_m22"],
			    m["m11"], m["m11"], m["m11"], m["m11"], m["m22"] );
  
  q["q11_3_q31"] = q_abcd( m["m11_3_m31"],
			   m["m11_3"], m["m11_2_m31"], m["m11_2_m31"], m["m11_2_m31"],
			   m["m11_2"], m["m11_2"], m["m11_m31"], m["m11_2"], m["m11_m31"], m["m11_m31"],
			   m["m11"], m["m11"], m["m11"], m["m31"] );

  q["q11_3_q32"] = q_abcd( m["m11_3_m32"],
			   m["m11_3"], m["m11_2_m32"], m["m11_2_m32"], m["m11_2_m32"],
			   m["m11_2"], m["m11_2"], m["m11_m32"], m["m11_2"], m["m11_m32"], m["m11_m32"],
			   m["m11"], m["m11"], m["m11"], m["m32"] );

  q["q11_3_q33"] = q_abcd( m["m11_3_m33"],
			   m["m11_3"], m["m11_2_m33"], m["m11_2_m33"], m["m11_2_m33"],
			   m["m11_2"], m["m11_2"], m["m11_m33"], m["m11_2"], m["m11_m33"], m["m11_m33"],
			   m["m11"], m["m11"], m["m11"], m["m33"] );
  
  q["q11_2_q22_q21"] = q_abcd( m["m11_2_m22_m21"],
			   m["m11_2_m22"], m["m11_2_m21"], m["m11_m21_m22"], m["m11_m21_m22"],
			   m["m11_2"], m["m11_m22"], m["m11_m21"], m["m11_m22"], m["m11_m21"], m["m21_m22"],
			       m["m11"], m["m11"], m["m22"], m["m21"] );


  q["q11_2_q21_2"] = q_abcd( m["m11_2_m21_2"],
			   m["m11_2_m21"], m["m11_2_m21"], m["m21_2_m11"], m["m21_2_m11"],
			   m["m11_2"], m["m11_m21"], m["m11_m21"], m["m11_m21"], m["m11_m21"], m["m21_2"],
			   m["m11"], m["m11"], m["m21"], m["m21"] );

  q["q11_2_q22_2"] = q_abcd( m["m11_2_m22_2"],
			   m["m11_2_m22"], m["m11_2_m22"], m["m22_2_m11"], m["m22_2_m11"],
			   m["m11_2"], m["m11_m22"], m["m11_m22"], m["m11_m22"], m["m11_m22"], m["m22_2"],
			   m["m11"], m["m11"], m["m22"], m["m22"] );

  
  q["q21_3"] = q_abc( m["m21_3"], m["m21_2"], m["m21_2"], m["m21_2"], m["m21"], m["m21"], m["m21"] );
  q["q22_3"] = q_abc( m["m22_3"], m["m22_2"], m["m22_2"], m["m22_2"], m["m22"], m["m22"], m["m22"] );

  q["q11_2_q41"] = q_abc( m["m11_2_m41"], m["m11_2"], m["m11_m41"], m["m11_m41"], m["m11"], m["m11"], m["m41"] );
  q["q11_2_q42"] = q_abc( m["m11_2_m42"], m["m11_2"], m["m11_m42"], m["m11_m42"], m["m11"], m["m11"], m["m42"] );
  q["q11_2_q43"] = q_abc( m["m11_2_m43"], m["m11_2"], m["m11_m43"], m["m11_m43"], m["m11"], m["m11"], m["m43"] );
  q["q11_2_q44"] = q_abc( m["m11_2_m44"], m["m11_2"], m["m11_m44"], m["m11_m44"], m["m11"], m["m11"], m["m44"] );

  q["q21_2_q22"] = q_abc( m["m21_2_m22"], m["m21_2"], m["m21_m22"], m["m21_m22"], m["m21"], m["m21"], m["m22"] );
  q["q22_2_q21"] = q_abc( m["m22_2_m21"], m["m22_2"], m["m21_m22"], m["m21_m22"], m["m22"], m["m22"], m["m21"] );


  q["q11_q21_q31"] = q_abc( m["m11_m21_m31"], m["m11_m21"], m["m11_m31"], m["m21_m31"], m["m11"], m["m21"], m["m31"] );
  q["q11_q21_q32"] = q_abc( m["m11_m21_m32"], m["m11_m21"], m["m11_m32"], m["m21_m32"], m["m11"], m["m21"], m["m32"] );
  q["q11_q21_q33"] = q_abc( m["m11_m21_m33"], m["m11_m21"], m["m11_m33"], m["m21_m33"], m["m11"], m["m21"], m["m33"] );
  q["q11_q22_q31"] = q_abc( m["m11_m22_m31"], m["m11_m22"], m["m11_m31"], m["m22_m31"], m["m11"], m["m22"], m["m31"] );
  q["q11_q22_q32"] = q_abc( m["m11_m22_m32"], m["m11_m22"], m["m11_m32"], m["m22_m32"], m["m11"], m["m22"], m["m32"] );
  q["q11_q22_q33"] = q_abc( m["m11_m22_m33"], m["m11_m22"], m["m11_m33"], m["m22_m33"], m["m11"], m["m22"], m["m33"] );

  
  q["q11_q51"] = q_ab( m["m11_m51"], m["m11"], m["m51"] );
  q["q11_q52"] = q_ab( m["m11_m52"], m["m11"], m["m52"] );
  q["q11_q53"] = q_ab( m["m11_m53"], m["m11"], m["m53"] );
  q["q11_q54"] = q_ab( m["m11_m54"], m["m11"], m["m54"] );
  q["q11_q55"] = q_ab( m["m11_m55"], m["m11"], m["m55"] );
  q["q21_q41"] = q_ab( m["m21_m41"], m["m21"], m["m41"] );
  q["q21_q42"] = q_ab( m["m21_m42"], m["m21"], m["m42"] );
  q["q21_q43"] = q_ab( m["m21_m43"], m["m21"], m["m43"] );
  q["q21_q44"] = q_ab( m["m21_m44"], m["m21"], m["m44"] );
  q["q22_q41"] = q_ab( m["m22_m41"], m["m22"], m["m41"] );
  q["q22_q42"] = q_ab( m["m22_m42"], m["m22"], m["m42"] );
  q["q22_q43"] = q_ab( m["m22_m43"], m["m22"], m["m43"] );
  q["q22_q44"] = q_ab( m["m22_m44"], m["m22"], m["m44"] );
  q["q31_2"]  = q_ab( m["m31_2"], m["m31"], m["m31"] );
  q["q31_q32"] = q_ab( m["m31_m32"], m["m32"], m["m32"] );
  q["q31_q33"] = q_ab( m["m31_m33"], m["m31"], m["m33"] );
  q["q32_2"] = q_ab( m["m32_2"], m["m32"], m["m32"] );
  q["q32_q33"] = q_ab( m["m32_m33"], m["m32"], m["m33"] );
  q["q33_2"] = q_ab( m["m33_2"], m["m33"], m["m33"] );
  
  q["q61"] = m["m61"];
  q["q62"] = m["m62"];
  q["q63"] = m["m63"];
  q["q64"] = m["m64"];
  q["q65"] = m["m65"];
  q["q66"] = m["m66"];

  //  cout << "m11_6 = " << m["m11_6"] << endl;

  q["q11_6"] = m["m11_6"] - 6*q["q11"]*q["q11_5"] - 15*q["q11_2"]*q["q11_4"] - 10*q["q11_3"]*q["q11_3"] - 15*q["q11"]*q["q11"]*q["q11_4"]
    - 60*q["q11"]*q["q11_2"]*q["q11_3"] - 15*pow(q["q11_2"],3) - 20*pow( q["q11"],3 )*q["q11_3"] -45*q["q11"]*q["q11"]*q["q11_2"]*q["q11_2"]
    -15*q["q11_2"]*pow(q["q11"],4) - pow(q["q11"],6);


}

//bivariate moments to cumulants
void CumulantProfileContainer::MomentsToCumulants()
{

  for(int i=0;i<300;i++)
    {
      double entries  = _profile["m11"]->GetBinEntries(i);
      if(entries < 1.0 ) continue;
      int cent = _profile["m11"]->GetBinCenter(i);
      
      InitCumulantMap(cent); // Fills all the m[] values
      
      double Q1 = q["q11"];
      
      double Q2 = q["q11_2"] + q["q21"] - q["q22"];

      double Q3 = q["q11_3"] + 3*q["q11_q21"] - 3*q["q11_q22"] +  q["q31"] - 3*q["q32"] + 2*q["q33"];

      double Q4 = q["q11_4"]
        + 6*q["q11_2_q21"] - 6*q["q11_2_q22"]
	+ 4*q["q11_q31"]
	+ 3*q["q21_2"] + 3*q["q22_2"]
       	-12*q["q11_q32"] + 8*q["q11_q33"]
	-6*q["q21_q22"]
	+ q["q41"] -7*q["q42"] + 12*q["q43"] - 6*q["q44"];

      double Q5 = q["q11_5"] + 10*q["q11_3_q21"] - 10*q["q11_3_q22"] + 10*q["q11_2_q31"] - 30*q["q11_2_q32"]
	+ 20*q["q11_2_q33"] + 15*q["q22_2_q11"] +15*q["q21_2_q11"] - 30*q["q11_q21_q22"]
	+5*q["q11_q41"] -35*q["q11_q42"] + 60*q["q11_q43"] -30*q["q11_q44"]
	+10*q["q21_q31"] - 30*q["q21_q32"] + 20*q["q21_q33"]
	-10*q["q22_q31"] +30*q["q22_q32"] -20*q["q22_q33"]
	+q["q51"] -15*q["q52"] +50*q["q53"] - 60*q["q54"] + 24*q["q55"];

      long double Q6 =  -120*pow(m["m11"],6)+360*pow(m["m11"],4)*m["m11_2"]-270*pow(m["m11"],2)*pow(m["m11_2"],2)
	+30*pow(m["m11_2"],3)-120*pow(m["m11"],3)*m["m11_3"]+120*m["m11"]*m["m11_2"]*m["m11_3"]
	-10*pow(m["m11_3"],2)+30*pow(m["m11"],2)*m["m11_4"]-15*m["m11_2"]*m["m11_4"]-6*m["m11"]*m["m11_5"]
	+m["m11_6"]
	+15*(24*pow(m["m11"],4)*m["m21"]-24*pow(m["m11"],3)*m["m11_m21"]-36*m["m21"]*pow(m["m11"],2)*m["m11_2"]
	     +24*m["m11"]*m["m11_m21"]*m["m11_2"]+6*m["m21"]*pow(m["m11_2"],2)+12*pow(m["m11"],2)*m["m11_2_m21"]
	     -6*m["m11_2"]*m["m11_2_m21"]+8*m["m21"]*m["m11"]*m["m11_3"]-4*m["m11_m21"]*m["m11_3"]
	     -4*m["m11"]*m["m11_3_m21"]-m["m21"]*m["m11_4"]+m["m11_4_m21"])
	-15*(24*pow(m["m11"],4)*m["m22"]-24*pow(m["m11"],3)*m["m11_m22"]-36*m["m22"]*pow(m["m11"],2)*m["m11_2"]
	     +24*m["m11"]*m["m11_m22"]*m["m11_2"]+6*m["m22"]*pow(m["m11_2"],2)+12*pow(m["m11"],2)*m["m11_2_m22"]
	     -6*m["m11_2"]*m["m11_2_m22"]+8*m["m22"]*m["m11"]*m["m11_3"]-4*m["m11_m22"]*m["m11_3"]
	     -4*m["m11"]*m["m11_3_m22"]-m["m22"]*m["m11_4"]+m["m11_4_m22"])
	+20*(-6*pow(m["m11"],3)*m["m31"]+6*m["m11"]*m["m11_2"]*m["m31"]-m["m11_3"]*m["m31"]
	     +6*pow(m["m11"],2)*m["m11_m31"]-3*m["m11_2"]*m["m11_m31"]-3*m["m11"]*m["m11_2_m31"]
	     +m["m11_3_m31"])
	-60*(-6*pow(m["m11"],3)*m["m32"]+6*m["m11"]*m["m11_2"]*m["m32"]-m["m11_3"]*m["m32"]
	     +6*pow(m["m11"],2)*m["m11_m32"]-3*m["m11_2"]*m["m11_m32"]-3*m["m11"]*m["m11_2_m32"]
	     +m["m11_3_m32"])
	+40*(-6*pow(m["m11"],3)*m["m33"]+6*m["m11"]*m["m11_2"]*m["m33"]-m["m11_3"]*m["m33"]
	     +6*pow(m["m11"],2)*m["m11_m33"]-3*m["m11_2"]*m["m11_m33"]-3*m["m11"]*m["m11_2_m33"]
	     +m["m11_3_m33"])
	-90*(-6*m["m21"]*m["m22"]*pow(m["m11"],2)+2*m["m21_m22"]*pow(m["m11"],2)+4*m["m22"]*m["m11"]*m["m11_m21"]
	     +4*m["m21"]*m["m11"]*m["m11_m22"]-2*m["m11_m21"]*m["m11_m22"]-2*m["m11"]*m["m11_m21_m22"]
	     +2*m["m21"]*m["m22"]*m["m11_2"]-m["m21_m22"]*m["m11_2"]-m["m22"]*m["m11_2_m21"]-m["m21"]*m["m11_2_m22"]
	     +m["m11_2_m22_m21"])
	+45*(-6*pow(m["m21"],2)*pow(m["m11"],2)+2*pow(m["m11"],2)*m["m21_2"]+8*m["m11"]*m["m21"]*m["m11_m21"]
	     -2*pow(m["m11_m21"],2)-2*m["m11"]*m["m21_2_m11"]+2*pow(m["m21"],2)*m["m11_2"]-m["m11_2"]*m["m21_2"]
	     -2*m["m21"]*m["m11_2_m21"]+m["m11_2_m21_2"])
	+45*(-6*pow(m["m22"],2)*pow(m["m11"],2)+2*pow(m["m11"],2)*m["m22_2"]+8*m["m11"]*m["m22"]*m["m11_m22"]
	     -2*pow(m["m11_m22"],2)-2*m["m11"]*m["m22_2_m11"]+2*pow(m["m22"],2)*m["m11_2"]-m["m11_2"]*m["m22_2"]
	     -2*m["m22"]*m["m11_2_m22"]+m["m11_2_m22_2"])
	+15*(2*pow(m["m21"],3)-3*m["m21"]*m["m21_2"]+m["m21_3"])
	-15*(2*pow(m["m22"],3)-3*m["m22"]*m["m22_2"]+m["m22_3"])
	+15* (m["m11_2_m41"]-m["m11_2"]*m["m41"]-2*m["m11"]*m["m11_m41"]+2*pow(m["m11"],2)*m["m41"])
	-105*(m["m11_2_m42"]-m["m11_2"]*m["m42"]-2*m["m11"]*m["m11_m42"]+2*pow(m["m11"],2)*m["m42"])
	+180*(m["m11_2_m43"]-m["m11_2"]*m["m43"]-2*m["m11"]*m["m11_m43"]+2*pow(m["m11"],2)*m["m43"])
	-90* (m["m11_2_m44"]-m["m11_2"]*m["m44"]-2*m["m11"]*m["m11_m44"]+2*pow(m["m11"],2)*m["m44"])
	-45* (m["m21_2_m22"]-m["m21_2"]*m["m22"]-2*m["m21"]*m["m21_m22"]+2*pow(m["m21"],2)*m["m22"])
	+45* (m["m22_2_m21"]-m["m22_2"]*m["m21"]-2*m["m22"]*m["m21_m22"]+2*pow(m["m22"],2)*m["m21"])
	+60* (m["m11_m21_m31"]-m["m11_m21"]*m["m31"]-m["m11_m31"]*m["m21"]-m["m21_m31"]*m["m11"]+2*m["m11"]*m["m21"]*m["m31"])
	-180*(m["m11_m21_m32"]-m["m11_m21"]*m["m32"]-m["m11_m32"]*m["m21"]-m["m21_m32"]*m["m11"]+2*m["m11"]*m["m21"]*m["m32"])
	+120*(m["m11_m21_m33"]-m["m11_m21"]*m["m33"]-m["m11_m33"]*m["m21"]-m["m21_m33"]*m["m11"]+2*m["m11"]*m["m21"]*m["m33"])
	-60* (m["m11_m22_m31"]-m["m11_m22"]*m["m31"]-m["m11_m31"]*m["m22"]-m["m22_m31"]*m["m11"]+2*m["m11"]*m["m22"]*m["m31"])
	+180*(m["m11_m22_m32"]-m["m11_m22"]*m["m32"]-m["m11_m32"]*m["m22"]-m["m22_m32"]*m["m11"]+2*m["m11"]*m["m22"]*m["m32"])
	-120*(m["m11_m22_m33"]-m["m11_m22"]*m["m33"]-m["m11_m33"]*m["m22"]-m["m22_m33"]*m["m11"]+2*m["m11"]*m["m22"]*m["m33"])
	+6*  (m["m11_m51"]-m["m11"]*m["m51"])
	-90* (m["m11_m52"]-m["m11"]*m["m52"])
	+300*(m["m11_m53"]-m["m11"]*m["m53"])
	-360*(m["m11_m54"]-m["m11"]*m["m54"])
	+144*(m["m11_m55"]-m["m11"]*m["m55"])
	+15* (m["m21_m41"]-m["m21"]*m["m41"])
	-105*(m["m21_m42"]-m["m21"]*m["m42"])
	+180*(m["m21_m43"]-m["m21"]*m["m43"])
	-90* (m["m21_m44"]-m["m21"]*m["m44"])
	-15* (m["m22_m41"]-m["m22"]*m["m41"])
	+105*(m["m22_m42"]-m["m22"]*m["m42"])
	-180*(m["m22_m43"]-m["m22"]*m["m43"])
	+90* (m["m22_m44"]-m["m22"]*m["m44"])
	+10* (m["m31_2"]  -  pow(m["m31"],2))
	-60* (m["m31_m32"]-m["m31"]*m["m32"])
	+40* (m["m31_m33"]-m["m31"]*m["m33"])
	+90* (m["m32_2"]  -  pow(m["m32"],2))
	-120*(m["m32_m33"]-m["m32"]*m["m33"])
	+40* (m["m33_2"]  -  pow(m["m33"],2))
	+    m["m61"] 
	-31* m["m62"] 
	+180*m["m63"] 
	-390*m["m64"] 
	+360*m["m65"]
	-120*m["m66"];
      
      double K1 = Q1;
      double K2 = Q2 - K1;
      double K3 = Q3 - 3*Q2 +2*Q1;
      double K4 = -6*Q1 + 11*Q2 -6*Q3 + Q4;
      double K5 = -10*Q4 + 35*Q3 -50*Q2 + 24*Q1 + Q5;
      double K6 = -120*Q1 + 274*Q2 - 225*Q3 + 85*Q4 - 15*Q5 + Q6;

      _gr["N"]->SetPoint(_gr["N"]->GetN(),cent, entries);
      _gr["C1"]->SetPoint(_gr["C1"]->GetN(), cent, Q1);
      _gr["C2"]->SetPoint(_gr["C2"]->GetN(), cent, Q2);
      _gr["C3"]->SetPoint(_gr["C3"]->GetN(), cent, Q3);
      _gr["C4"]->SetPoint(_gr["C4"]->GetN(), cent, Q4);
      _gr["C5"]->SetPoint(_gr["C5"]->GetN(), cent, Q5);
      _gr["C6"]->SetPoint(_gr["C6"]->GetN(), cent, Q6);

      _gr["K1"]->SetPoint(_gr["K1"]->GetN(), cent, K1);
      _gr["K2"]->SetPoint(_gr["K2"]->GetN(), cent, K2);
      _gr["K3"]->SetPoint(_gr["K3"]->GetN(), cent, K3);
      _gr["K4"]->SetPoint(_gr["K4"]->GetN(), cent, K4);
      _gr["K5"]->SetPoint(_gr["K5"]->GetN(), cent, K5);
      _gr["K6"]->SetPoint(_gr["K6"]->GetN(), cent, K6);

    }
  
}

//Calls ReBin( ... ) for all cumulants and cumulant ratios
void CumulantProfileContainer::ReBinAllGraphs(std::vector<int> binEdges, std::vector<double> binLabels)
{

  ReBin(binEdges,binLabels,_gr["N"],_gr["C1"],_gr["C1_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["C2"],_gr["C2_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["C3"],_gr["C3_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["C4"],_gr["C4_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["C5"],_gr["C5_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["C6"],_gr["C6_cbwc"]);

  ReBin(binEdges,binLabels,_gr["N"],_gr["K1"],_gr["K1_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["K2"],_gr["K2_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["K3"],_gr["K3_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["K4"],_gr["K4_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["K5"],_gr["K5_cbwc"]);
  ReBin(binEdges,binLabels,_gr["N"],_gr["K6"],_gr["K6_cbwc"]);

  RatioGraphs( _gr["C2"], _gr["C1"], _gr["C2_C1"] );
  RatioGraphs( _gr["C3"], _gr["C2"], _gr["C3_C2"] );
  RatioGraphs( _gr["C4"], _gr["C2"], _gr["C4_C2"] );
  RatioGraphs( _gr["C5"], _gr["C1"], _gr["C5_C1"] );
  RatioGraphs( _gr["C6"], _gr["C2"], _gr["C6_C2"] );

  RatioGraphs( _gr["K2"], _gr["K1"], _gr["K2_K1"] );
  RatioGraphs( _gr["K3"], _gr["K1"], _gr["K3_K1"] );
  RatioGraphs( _gr["K4"], _gr["K1"], _gr["K4_K1"] );
  RatioGraphs( _gr["K5"], _gr["K1"], _gr["K5_K1"] );
  RatioGraphs( _gr["K6"], _gr["K1"], _gr["K6_K1"] );

  RatioGraphs( _gr["C2_cbwc"], _gr["C1_cbwc"], _gr["C2_C1_cbwc"] );
  RatioGraphs( _gr["C3_cbwc"], _gr["C2_cbwc"], _gr["C3_C2_cbwc"] );
  RatioGraphs( _gr["C4_cbwc"], _gr["C2_cbwc"], _gr["C4_C2_cbwc"] );
  RatioGraphs( _gr["C5_cbwc"], _gr["C1_cbwc"], _gr["C5_C1_cbwc"] );
  RatioGraphs( _gr["C6_cbwc"], _gr["C2_cbwc"], _gr["C6_C2_cbwc"] );

  RatioGraphs( _gr["K2_cbwc"], _gr["K1_cbwc"], _gr["K2_K1_cbwc"] );
  RatioGraphs( _gr["K3_cbwc"], _gr["K1_cbwc"], _gr["K3_K1_cbwc"] );
  RatioGraphs( _gr["K4_cbwc"], _gr["K1_cbwc"], _gr["K4_K1_cbwc"] );
  RatioGraphs( _gr["K5_cbwc"], _gr["K1_cbwc"], _gr["K5_K1_cbwc"] );
  RatioGraphs( _gr["K6_cbwc"], _gr["K1_cbwc"], _gr["K6_K1_cbwc"] );
  
}

//Perform CBWC
void CumulantProfileContainer::ReBin(std::vector<int> &binEdges,std::vector<double> &binLabels,TGraphErrors * Ngr, TGraphErrors * Qgr,TGraphErrors * Qgr_cbwc)
{

  int nbins = binEdges.size() - 1;
  if (binLabels.size() != nbins)
    {
      std::cout << "binLabel array size (size() -1) doesnt match binEdge size" << std::endl;
    }

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> w;
  
  for (int i=0; i<Ngr->GetN(); i++)
    {
      x.push_back(Ngr->GetX()[i]);
      y.push_back(Qgr->GetY()[i]);
      w.push_back(Ngr->GetY()[i]);
    }
  
  for (int i=0;i<nbins;i++)
    {

      double weighted_sum = 0;
      double weight = 0;
      
      for (int j=0;j<x.size();j++)
	{
	  if ( x[j] >= binEdges[i] && x[j] < binEdges[i+1] )
		{
	      weighted_sum += y[j]*w[j];
	      weight += w[j];
	    }
	}

      double val = weighted_sum/weight;
      Qgr_cbwc->SetPoint(Qgr_cbwc->GetN(), binLabels[i], val );
      
    }
  
}

//Take the ratio of two graphs
void CumulantProfileContainer::RatioGraphs(TGraphErrors * num_gr, TGraphErrors * denom_gr, TGraphErrors * ratioGr)
{
  if ( num_gr->GetN() != denom_gr->GetN())
    {
      std::cerr << "Graphs have different number of points";
      exit(1);
    }
  
  for (int iPoint=0; iPoint<num_gr->GetN();iPoint++)
    {
      if (fabs(num_gr->GetX()[iPoint] - denom_gr->GetX()[iPoint]) > 0) continue;
      double ratio =num_gr->GetY()[iPoint]/denom_gr->GetY()[iPoint];
      ratioGr->SetPoint(ratioGr->GetN(),num_gr->GetX()[iPoint],ratio);
    }
}
