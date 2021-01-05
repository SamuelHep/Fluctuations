#include <iostream>
#include <vector>
#include <utility>
#include <map>

#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "CumulantProfileContainer.h"

ClassImp(CumulantProfileContainer)

using namespace std;

CumulantProfileContainer::CumulantProfileContainer()
{

  for (int i=0;i<23;i++)
    {
      TString prof_name = TString::Format("Profile_%s",IndexToName(i).Data());
      TProfile * temp = new TProfile(prof_name.Data(),prof_name.Data(),300,-0.5,299.5);
      _prof_vec.push_back(temp);
    }

  _gr_N = new TGraphErrors();  _gr_N->SetName("Ngr");
  _gr_Q1 = new TGraphErrors(); _gr_Q1->SetName("Q1gr");
  _gr_Q2 = new TGraphErrors(); _gr_Q2->SetName("Q2gr");
  _gr_Q3 = new TGraphErrors(); _gr_Q3->SetName("Q3gr");
  _gr_Q4 = new TGraphErrors(); _gr_Q4->SetName("Q4gr");
  _gr_binned_Q1 = NULL;
  _gr_binned_Q2 = NULL;
  _gr_binned_Q3 = NULL;
  _gr_binned_Q4 = NULL;
  _gr_Q2_Q1 = NULL;
  _gr_Q3_Q2 = NULL;
  _gr_Q4_Q2 = NULL;

  _gr_binned_K1 = NULL;
  _gr_binned_K2 = NULL;
  _gr_binned_K3 = NULL;
  _gr_binned_K4 = NULL;
  _gr_K2_K1 = NULL;
  _gr_K3_K1 = NULL;
  _gr_K4_K1 = NULL;
  
  _gr_vec.push_back(_gr_N);
  _gr_vec.push_back(_gr_Q1);
  _gr_vec.push_back(_gr_Q2);
  _gr_vec.push_back(_gr_Q3);
  _gr_vec.push_back(_gr_Q4);

  _qGraphsFilled = false;

}

CumulantProfileContainer::CumulantProfileContainer(int iName)
{
  for (int i=0;i<23;i++)
    {
      TString prof_name = TString::Format("Profile_%s_%i",IndexToName(i).Data(),iName);
      TProfile * temp = new TProfile(prof_name.Data(),prof_name.Data(),300,-0.5,299.5);
      _prof_vec.push_back(temp);
    }

  _gr_N = new TGraphErrors();  _gr_N->SetName(TString::Format("Ngr_%i",iName));
  _gr_Q1 = new TGraphErrors(); _gr_Q1->SetName(TString::Format("Q1gr_%i",iName));
  _gr_Q2 = new TGraphErrors(); _gr_Q2->SetName(TString::Format("Q2gr_%i",iName));
  _gr_Q3 = new TGraphErrors(); _gr_Q3->SetName(TString::Format("Q3gr_%i",iName));
  _gr_Q4 = new TGraphErrors(); _gr_Q4->SetName(TString::Format("Q4gr_%i",iName));
  _gr_binned_Q1 = NULL;
  _gr_binned_Q2 = NULL;
  _gr_binned_Q3 = NULL;
  _gr_binned_Q4 = NULL;

  _gr_binned_K1 = NULL;
  _gr_binned_K2 = NULL;
  _gr_binned_K3 = NULL;
  _gr_binned_K4 = NULL;
  _gr_K2_K1 = NULL;
  _gr_K3_K1 = NULL;
  _gr_K4_K1 = NULL;

  _gr_vec.push_back(_gr_N);
  _gr_vec.push_back(_gr_Q1);
  _gr_vec.push_back(_gr_Q2);
  _gr_vec.push_back(_gr_Q3);
  _gr_vec.push_back(_gr_Q4);

  _qGraphsFilled = false;
  
}

CumulantProfileContainer::CumulantProfileContainer(int iName,TString label,TH1F * multHist,LongDouble_t Cumulants[300][5],int MaxMult)
{

  if ( iName == kPrimary )
    {
      _gr_N = new TGraphErrors();  _gr_N->SetName(TString::Format("Ngr_%s",label.Data()));
      _gr_Q1 = new TGraphErrors(); _gr_Q1->SetName(TString::Format("Q1gr_%s",label.Data()));
      _gr_Q2 = new TGraphErrors(); _gr_Q2->SetName(TString::Format("Q2gr_%s",label.Data()));
      _gr_Q3 = new TGraphErrors(); _gr_Q3->SetName(TString::Format("Q3gr_%s",label.Data()));
      _gr_Q4 = new TGraphErrors(); _gr_Q4->SetName(TString::Format("Q4gr_%s",label.Data()));
    }
  else 
    {
      _gr_N = new TGraphErrors();  _gr_N->SetName(TString::Format("Ngr_%i_%s",iName,label.Data()));
      _gr_Q1 = new TGraphErrors(); _gr_Q1->SetName(TString::Format("Q1gr_%i_%s",iName,label.Data()));
      _gr_Q2 = new TGraphErrors(); _gr_Q2->SetName(TString::Format("Q2gr_%i_%s",iName,label.Data()));
      _gr_Q3 = new TGraphErrors(); _gr_Q3->SetName(TString::Format("Q3gr_%i_%s",iName,label.Data()));
      _gr_Q4 = new TGraphErrors(); _gr_Q4->SetName(TString::Format("Q4gr_%i_%s",iName,label.Data()));
    }

  _gr_binned_Q1 = NULL;
  _gr_binned_Q2 = NULL;
  _gr_binned_Q3 = NULL;
  _gr_binned_Q4 = NULL;

  _gr_binned_K1 = NULL;
  _gr_binned_K2 = NULL;
  _gr_binned_K3 = NULL;
  _gr_binned_K4 = NULL;
  _gr_K2_K1 = NULL;
  _gr_K3_K1 = NULL;
  _gr_K4_K1 = NULL;

  for (int iMult=0;iMult<MaxMult;iMult++)
    {
      _gr_N->SetPoint(_gr_N->GetN(),iMult,multHist->GetBinContent(multHist->FindBin(iMult)));
      _gr_Q1->SetPoint(_gr_Q1->GetN(),iMult,Cumulants[iMult][1]);
      _gr_Q2->SetPoint(_gr_Q2->GetN(),iMult,Cumulants[iMult][2]);
      _gr_Q3->SetPoint(_gr_Q3->GetN(),iMult,Cumulants[iMult][3]);
      _gr_Q4->SetPoint(_gr_Q4->GetN(),iMult,Cumulants[iMult][4]);
    }

  _gr_vec.push_back(_gr_N);
  _gr_vec.push_back(_gr_Q1);
  _gr_vec.push_back(_gr_Q2);
  _gr_vec.push_back(_gr_Q3);
  _gr_vec.push_back(_gr_Q4);

  _qGraphsFilled = true;
  
}


CumulantProfileContainer::CumulantProfileContainer(TFile *f,int iName)
{

  if ( iName == kPrimary )
    {
      for (int i=0;i<23;i++)
	{
	  TString prof_name = TString::Format("Profile_%s",IndexToName(i).Data());
	  TProfile * temp = (TProfile*) f->Get(prof_name);
	  _prof_vec.push_back(temp);
	}

      _gr_N = new TGraphErrors();  _gr_N->SetName("Ngr");
      _gr_Q1 = new TGraphErrors(); _gr_Q1->SetName("Q1gr");
      _gr_Q2 = new TGraphErrors(); _gr_Q2->SetName("Q2gr");
      _gr_Q3 = new TGraphErrors(); _gr_Q3->SetName("Q3gr");
      _gr_Q4 = new TGraphErrors(); _gr_Q4->SetName("Q4gr");
      
    }
  else{
    for (int i=0;i<23;i++)
      {
	TString prof_name = TString::Format("Profile_%s_%i",IndexToName(i).Data(),iName);
	TProfile * temp = (TProfile*) f->Get(prof_name);
	_prof_vec.push_back(temp);
      }
  
    _gr_N = new TGraphErrors();  _gr_N->SetName(TString::Format("Ngr_%i",iName));
    _gr_Q1 = new TGraphErrors(); _gr_Q1->SetName(TString::Format("Q1gr_%i",iName));
    _gr_Q2 = new TGraphErrors(); _gr_Q2->SetName(TString::Format("Q2gr_%i",iName));
    _gr_Q3 = new TGraphErrors(); _gr_Q3->SetName(TString::Format("Q3gr_%i",iName));
    _gr_Q4 = new TGraphErrors(); _gr_Q4->SetName(TString::Format("Q4gr_%i",iName));
  }
  
  _gr_binned_Q1 = NULL;
  _gr_binned_Q2 = NULL;
  _gr_binned_Q3 = NULL;
  _gr_binned_Q4 = NULL;

  _gr_binned_K1 = NULL;
  _gr_binned_K2 = NULL;
  _gr_binned_K3 = NULL;
  _gr_binned_K4 = NULL;
  _gr_K2_K1 = NULL;
  _gr_K3_K1 = NULL;
  _gr_K4_K1 = NULL;

  _gr_vec.push_back(_gr_N);
  _gr_vec.push_back(_gr_Q1);
  _gr_vec.push_back(_gr_Q2);
  _gr_vec.push_back(_gr_Q3);
  _gr_vec.push_back(_gr_Q4);

  _qGraphsFilled = false;

}

CumulantProfileContainer::~CumulantProfileContainer()
{

  for (auto &v : _gr_vec )
    {
      if (!v) continue;
      v->Delete();
    }

  _gr_vec.clear();

  for (auto &v : _prof_vec )
    {
      if (!v) continue;
      v->Delete();
    }

  _prof_vec.clear();
}


TGraphErrors * CumulantProfileContainer::GetNGraph(int n)
{
  if ( !_qGraphsFilled ) FactorialCumulantsToNormalCumulants();
  return _gr_vec[n];
}


void CumulantProfileContainer::FillProfile(int cent,std::vector<std::vector<double>> m_r_s)
{

  FillFactorialCumulantQuants(m_r_s);
  for (unsigned int i=0;i<23;i++)
    {
      _prof_vec[i]->Fill(cent,_m[i]);
    }
  
}

void CumulantProfileContainer::FactorialCumulantsToNormalCumulants()
{

  for(int i=0;i<300;i++)
    {
     double entries  = _prof_vec[0]->GetBinEntries(i);
      if(entries < 1.0 ) continue;
      int cent = _prof_vec[0]->GetBinCenter(i);

      std::map<TString,double> q = InitMap(cent);
            
      double Q1 = q["q11"];
      
      double Q2 = q["q11q11"] + q["q21"] - q["q22"];
      
      double Q3 = q["q11q11q11"] + 3*q["q11q21"] - 3*q["q11q22"] +  q["q31"] - 3*q["q32"] + 2*q["q33"];

      double Q4 = q["q11q11q11q11"]
        + 6*q["q11q11q21"] - 6*q["q11q11q22"]
	+ 4*q["q11q31"]
	+ 3*q["q21q21"] + 3*q["q22q22"]
       	-12*q["q11q32"] + 8*q["q11q33"]
	-6*q["q21q22"]
	+ q["q41"] -7*q["q42"] + 12*q["q43"] - 6*q["q44"];


      _gr_N->SetPoint(_gr_N->GetN(),cent, entries);
      _gr_Q1->SetPoint(_gr_Q1->GetN(),cent, Q1);
      _gr_Q2->SetPoint(_gr_Q2->GetN(),cent, Q2);
      _gr_Q3->SetPoint(_gr_Q3->GetN(),cent, Q3);
      _gr_Q4->SetPoint(_gr_Q4->GetN(),cent, Q4);

    }

  _qGraphsFilled = true;
  
}

void CumulantProfileContainer::ReBinAllGraphs(std::vector<int> binEdges, std::vector<double> binLabels)
{

  _gr_binned_Q1 = ReBin(binEdges,binLabels,_gr_N,_gr_Q1);
  _gr_binned_Q2 = ReBin(binEdges,binLabels,_gr_N,_gr_Q2);
  _gr_binned_Q3 = ReBin(binEdges,binLabels,_gr_N,_gr_Q3);
  _gr_binned_Q4 = ReBin(binEdges,binLabels,_gr_N,_gr_Q4);

  _gr_Q2_Q1 = RatioGraphs(_gr_binned_Q2,_gr_binned_Q1);
  _gr_Q3_Q2 = RatioGraphs(_gr_binned_Q3,_gr_binned_Q2);
  _gr_Q4_Q2 = RatioGraphs(_gr_binned_Q4,_gr_binned_Q2);
  
  _gr_binned_K1 = CalcK( 1 );
  _gr_binned_K2 = CalcK( 2 );
  _gr_binned_K3 = CalcK( 3 );
  _gr_binned_K4 = CalcK( 4 );

  _gr_K2_K1 = RatioGraphs(_gr_binned_K2,_gr_binned_K1);
  _gr_K3_K1 = RatioGraphs(_gr_binned_K3,_gr_binned_K1);
  _gr_K4_K1 = RatioGraphs(_gr_binned_K4,_gr_binned_K1);

  _gr_vec.push_back(_gr_binned_Q1);
  _gr_vec.push_back(_gr_binned_Q2);
  _gr_vec.push_back(_gr_binned_Q3);
  _gr_vec.push_back(_gr_binned_Q4);

  _gr_vec.push_back(_gr_Q2_Q1);
  _gr_vec.push_back(_gr_Q3_Q2);
  _gr_vec.push_back(_gr_Q4_Q2);

  _gr_vec.push_back( _gr_binned_K1);
  _gr_vec.push_back( _gr_binned_K2);
  _gr_vec.push_back( _gr_binned_K3);
  _gr_vec.push_back( _gr_binned_K4);

  _gr_vec.push_back( _gr_K2_K1);
  _gr_vec.push_back( _gr_K3_K1);
  _gr_vec.push_back( _gr_K4_K1);

  
}


TGraphErrors * CumulantProfileContainer::CalcK(int i)
{

  if ( !_gr_binned_Q1 || !_gr_binned_Q2 || !_gr_binned_Q3 || !_gr_binned_Q4 )
    {
      std::cout << "_gr_binned_Q not filled" << std::endl;
      return NULL;
    }

  if ( i==1 ) return CalcK1();
  if ( i==2 ) return CalcK2();
  if ( i==3 ) return CalcK3();
  if ( i==4 ) return CalcK4();

  return NULL;
}

TGraphErrors * CumulantProfileContainer::CalcK1()
{
  //K1 = C1
  TGraphErrors * gr = new TGraphErrors();

  TString name = _gr_binned_Q1->GetName();
  name.ReplaceAll("Q","K");
  gr->SetName(name);

  if ( _gr_binned_Q1->GetN() == 0 )
    {
      cout << "Q1->GetN() == 0" << endl;
      return NULL;
    }

  for ( int i=0;i<_gr_binned_Q1->GetN(); i++)
    {
      double K1 = _gr_binned_Q1->GetY()[i];
      double x = _gr_binned_Q1->GetX()[i];
      double xerr = _gr_binned_Q1->GetX()[i];
      double yerr = _gr_binned_Q1->GetY()[i];

      gr->SetPoint(gr->GetN(),x,K1);
      gr->SetPointError(gr->GetN()-1,xerr,yerr);
    }

  return gr;
}

TGraphErrors * CumulantProfileContainer::CalcK2()
{

  //K2  = Q2 - Q1

  TGraphErrors * gr = new TGraphErrors();

  TString name = _gr_binned_Q2->GetName();
  name.ReplaceAll("Q","K");
  gr->SetName(name);

  if ( _gr_binned_Q1->GetN() != _gr_binned_Q2->GetN() )
    {
      cout << "Q1->GetN() != Q2->GetN()" << endl;
      return NULL;
    }

  for ( int i=0;i<_gr_binned_Q1->GetN(); i++)
    {
      double K1 = _gr_binned_Q1->GetY()[i];
      double Q2 = _gr_binned_Q2->GetY()[i];

      double x = _gr_binned_Q1->GetX()[i];
      double xerr = _gr_binned_Q1->GetX()[i];
      double yerr = _gr_binned_Q1->GetY()[i];

      double K2 = Q2 - K1;

      gr->SetPoint(gr->GetN(),x,K2);
      gr->SetPointError(gr->GetN()-1,xerr,yerr);
    }

  return gr;
}

TGraphErrors * CumulantProfileContainer::CalcK3()
{

  //K3  = Q3 - K1 - 3*K2

  TGraphErrors * gr = new TGraphErrors();

  TString name = _gr_binned_Q3->GetName();
  name.ReplaceAll("Q","K");
  gr->SetName(name);


  if ( _gr_binned_Q1->GetN() != _gr_binned_Q2->GetN() )
    {
      cout << "Q1->GetN() != Q2->GetN()" << endl;
      return NULL;
    }
  if ( _gr_binned_Q1->GetN() != _gr_binned_Q3->GetN() )
    {
      cout << "Q1->GetN() != Q3->GetN()" << endl;
      return NULL;
    }
  
  for ( int i=0;i<_gr_binned_Q3->GetN(); i++)
    {
      double Q1 = _gr_binned_Q1->GetY()[i];
      double Q2 = _gr_binned_Q2->GetY()[i];
      double Q3 = _gr_binned_Q3->GetY()[i];
      double x = _gr_binned_Q1->GetX()[i];
      double xerr = _gr_binned_Q1->GetX()[i];
      double yerr = _gr_binned_Q1->GetY()[i];

      double K3 = Q3 - 3*Q2 + 2*Q1;

      gr->SetPoint(gr->GetN(),x,K3);
      gr->SetPointError(gr->GetN()-1,xerr,yerr);
    }
  
  return gr;
}

TGraphErrors * CumulantProfileContainer::CalcK4()
{

  //K4 = Q4 - K1 - 7*K2 - 6*K3

  TGraphErrors * gr = new TGraphErrors();

  TString name = _gr_binned_Q4->GetName();
  name.ReplaceAll("Q","K");
  gr->SetName(name);

  if ( _gr_binned_Q1->GetN() != _gr_binned_Q2->GetN() )
    {
      cout << "Q1->GetN() != Q2->GetN()" << endl;
      return NULL;
    }
  if ( _gr_binned_Q1->GetN() != _gr_binned_Q3->GetN() )
    {
      cout << "Q1->GetN() != Q3->GetN()" << endl;
      return NULL;
    }
  if ( _gr_binned_Q1->GetN() != _gr_binned_Q4->GetN() )
    {
      cout << "Q1->GetN() != Q4->GetN()" << endl;
      return NULL;
    }

  for ( int i=0;i<_gr_binned_Q4->GetN(); i++)
    {
      double Q1 = _gr_binned_Q1->GetY()[i];
      double Q2 = _gr_binned_Q2->GetY()[i];
      double Q3 = _gr_binned_Q3->GetY()[i];
      double Q4 = _gr_binned_Q4->GetY()[i];
      double x = _gr_binned_Q1->GetX()[i];
      double xerr = _gr_binned_Q1->GetX()[i];
      double yerr = _gr_binned_Q1->GetY()[i];

      //      double K2 = Q2 - Q1;
      //      double K3 = Q3 - Q1 - -3*K2;
      //      double K4 = Q4 - Q1 - 7*K2 - 6*K3;

      double K4 = -6*Q1 + 11*Q2 -6*Q3 + Q4;

      gr->SetPoint(gr->GetN(),x,K4);
      gr->SetPointError(gr->GetN()-1,xerr,yerr);
    }

  return gr;
}

TGraphErrors * CumulantProfileContainer::RatioGraphs(TGraphErrors * num_gr, TGraphErrors * denom_gr)
{
  if ( num_gr->GetN() != denom_gr->GetN())
    {
      std::cerr << "Graphs have different number of points";
      exit(1);
    }

  TGraphErrors * ratioGr = new TGraphErrors();
  TString name = TString::Format("ratio_%s_%s",num_gr->GetName(),denom_gr->GetName());
  ratioGr->SetName(name);
  
  for (int iPoint=0; iPoint<num_gr->GetN();iPoint++)
    {
      if (fabs(num_gr->GetX()[iPoint] - denom_gr->GetX()[iPoint]) > 0) continue;
      double ratio =num_gr->GetY()[iPoint]/denom_gr->GetY()[iPoint];
      ratioGr->SetPoint(ratioGr->GetN(),num_gr->GetX()[iPoint],ratio);
    }

  return ratioGr;
}

//Provide both low and high edge... binEdges.size() == nbins + 1 
TGraphErrors * CumulantProfileContainer::ReBin(std::vector<int> &binEdges,std::vector<double> &binLabels,TGraphErrors * Ngr, TGraphErrors * Qgr)
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
  
  TGraphErrors * BinnedGr = new TGraphErrors();
  
  //Make sure if has a new name to avoid root name problems
  TString name = Qgr->GetName();
  TString binnedname = TString::Format("%s_binnned", name.Data()); 
  BinnedGr->SetName(binnedname);
  
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
      BinnedGr->SetPoint(BinnedGr->GetN(), binLabels[i], val );
      
    }
  
  return BinnedGr;
  
}

double CumulantProfileContainer::Avg_M(int cent,int index)
{
  return (double) _prof_vec[index]->GetBinContent(_prof_vec[index]->FindBin(cent)); // index 0 -> < q(1,1) > // index 1 -> < q(1,1)^2 >
}


TString CumulantProfileContainer::IndexToName(int i)
{

  TString name;
  
  if( i==0 ) name = "q11";
  if( i==1 ) name = "q11q11";
  if( i==2 ) name = "q21";
  if( i==3 ) name = "q22";
  if( i==4 ) name = "q11q11q11";
  if( i==5 ) name = "q11q21";
  if( i==6 ) name = "q11q22";
  if( i==7 ) name = "q31";
  if( i==8 ) name = "q32";
  if( i==9 ) name = "q33";
  if( i==10 ) name = "q11q11q11q11";
  if( i==11 ) name = "q11q11q21";
  if( i==12 ) name = "q11q11q22";
  if( i==13 ) name = "q11q31";
  if( i==14 ) name = "q21q21";
  if( i==15 ) name = "q22q22";
  if( i==16 ) name = "q11q32";
  if( i==17 ) name = "q11q33";
  if( i==18 ) name = "q21q22";
  if( i==19 ) name = "q41";
  if( i==20 ) name = "q42";
  if( i==21 ) name = "q43";
  if( i==22 ) name = "q44";
  if(i<0||i>22) name= "bad_index";

  return name;
}

void CumulantProfileContainer::FillFactorialCumulantQuants(std::vector<std::vector<double>> & m_r_s)
{

  _m[0] = m_r_s[1][1];
  _m[1] = m_r_s[1][1]*m_r_s[1][1]; //  q(1,1)^2 
  _m[2] = m_r_s[2][1]; //  q(2,1)^2 
  _m[3] = m_r_s[2][2]; 
  _m[4] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1];
  _m[5] = m_r_s[1][1]*m_r_s[2][1];
  _m[6] = m_r_s[1][1]*m_r_s[2][2];
  _m[7] = m_r_s[3][1];
  _m[8] = m_r_s[3][2];
  _m[9] = m_r_s[3][3];
  _m[10] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1]*m_r_s[1][1];
  _m[11] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][1];
  _m[12] = m_r_s[1][1]*m_r_s[1][1]*m_r_s[2][2];
  _m[13] = m_r_s[1][1]*m_r_s[3][1];
  _m[14] = m_r_s[2][1]*m_r_s[2][1];
  _m[15] = m_r_s[2][2]*m_r_s[2][2];
  _m[16] = m_r_s[1][1]*m_r_s[3][2];
  _m[17] = m_r_s[1][1]*m_r_s[3][3];
  _m[18] = m_r_s[2][1]*m_r_s[2][2];
  _m[19] = m_r_s[4][1];
  _m[20] = m_r_s[4][2];
  _m[21] = m_r_s[4][3];
  _m[22] = m_r_s[4][4];

}

std::map<TString,double> CumulantProfileContainer::InitMap(int cent)
{

  std::map<TString,double> q;
  
  long double m11 = Avg_M(cent,0);
  long double m11m11 = Avg_M(cent,1);
  long double m21 = Avg_M(cent,2);
  long double m22 = Avg_M(cent,3);
  long double m11m11m11 = Avg_M(cent,4);
  long double m11m21 = Avg_M(cent,5);
  long double m11m22 = Avg_M(cent,6);
  long double m31 = Avg_M(cent,7);
  long double m32 = Avg_M(cent,8);
  long double m33 = Avg_M(cent,9);
  long double m11m11m11m11 = Avg_M(cent,10);
  long double m11m11m21 = Avg_M(cent,11);
  long double m11m11m22 = Avg_M(cent,12);
  long double m11m31 = Avg_M(cent,13);
  long double m21m21 = Avg_M(cent,14);
  long double m22m22 = Avg_M(cent,15);
  long double m11m32 = Avg_M(cent,16);
  long double m11m33 = Avg_M(cent,17);
  long double m21m22 = Avg_M(cent,18);
  long double m41 = Avg_M(cent,19);
  long double m42 = Avg_M(cent,20);
  long double m43 = Avg_M(cent,21);
  long double m44 = Avg_M(cent,22);

  //Convert the moments to cumulants 
  
  long double q11         = m11;
  long double q11q11      = m11m11 - m11*m11;
  long double q21         = m21;
  long double q22         = m22; 
  long double q11q11q11   = m11m11m11 - 3*m11m11*m11 + 2*m11*m11*m11;
  long double q11q21      = m11m21 - m11*m21;  
  long double q11q22      = m11m22 - m11*m22;
  long double q31         = m31;
  long double q32         = m32;
  long double q33         = m33;
  long double q11q11q11q11= m11m11m11m11 - 4*m11m11m11*m11 -3*m11m11*m11m11 + //q(1,1)^4
                     12*m11m11*m11*m11 -6*m11*m11*m11*m11;
  long double q11q11q21   = m11m11m21 -  q11q11*m21 - 2*q11q21*m11 - m11*m11*m21;
  long double q11q11q22   = m11m11m22 -  q11q11*m22 - 2*q11q22*m11 - m11*m11*m22;
  long double q11q31      = m11m31 - m11*m31;
  long double q21q21      = m21m21 - m21*m21;
  long double q22q22      = m22m22 - m22*m22;
  long double q11q32      = m11m32 - m11*m32;
  long double q11q33      = m11m33 - m11*m33;
  long double q21q22      = m21m22 - m21*m22;
  long double q41         = m41;
  long double q42         = m42;
  long double q43         = m43;
  long double q44         = m44; // q(4,4)
  
  q["q11"] = q11;
  q["q11q11"] = q11q11;
  q["q21"] =q21;
  q["q22"] = q22;
  q["q11q11q11"] = q11q11q11;
  q["q11q21"] = q11q21;
  q["q11q22"] = q11q22;
  q["q31"] = q31;
  q["q32"] = q32;
  q["q33"] = q33;
  q["q11q11q11q11"] = q11q11q11q11;
  q["q11q11q21"] = q11q11q21;
  q["q11q11q22"] = q11q11q22;
  q["q11q31"] = q11q31;
  q["q21q21"] = q21q21;
  q["q22q22"] = q22q22;
  q["q11q32"] = q11q32;
  q["q11q33"] = q11q33;
  q["q21q22"] = q21q22;
  q["q41"] = q41;
  q["q42"] = q42;
  q["q43"] = q43;
  q["q44"] = q44;
  
  return q;
   
}
