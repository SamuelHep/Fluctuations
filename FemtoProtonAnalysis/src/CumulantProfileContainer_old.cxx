#include <iostream>
#include <vector>
#include <utility>
#include <map>

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

      cout << "Q3=" << Q3 << endl;
      cout << " q11q11q11=" << q["q11q11q11"] << endl;
      cout << " q11q21=" << q["q11q21"] << endl;
      cout << " q11q22=" << q["q11q22"] << endl;
      cout << " q11q31=" << q["q11q31"] << endl;
      cout << " q31=" << q["q31"] << endl;
      cout << " q32=" << q["q32"] << endl;
      cout << " q33=" << q["q33"] << endl;
      
      double Q4 = q["q11q11q11q11"]
        + 6*q["q11q21"] - 6*q["q11q22"]
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
  
  _gr_vec.push_back(_gr_binned_Q1);
  _gr_vec.push_back(_gr_binned_Q2);
  _gr_vec.push_back(_gr_binned_Q3);
  _gr_vec.push_back(_gr_binned_Q4);

  _gr_vec.push_back(_gr_Q2_Q1);
  _gr_vec.push_back(_gr_Q3_Q2);
  _gr_vec.push_back(_gr_Q4_Q2);
  
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

  std::cout << "DEBUG::0.1" << std::endl;
  
  int nbins = binEdges.size() - 1;
  if (binLabels.size() != nbins)
    {
      std::cout << "binLabel array size (size() -1) doesnt match binEdge size" << std::endl;
    }

  std::cout << "DEBUG::0.2" << std::endl;
  
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> w;
  
  for (int i=0; i<Ngr->GetN(); i++)
    {
      x.push_back(Ngr->GetX()[i]);
      y.push_back(Qgr->GetY()[i]);
      w.push_back(Ngr->GetY()[i]);
    }

  std::cout << "DEBUG::0.3" << std::endl;
  
  TGraphErrors * BinnedGr = new TGraphErrors();

  std::cout << "DEBUG::0.4" << std::endl;
  
  //Make sure if has a new name to avoid root name problems
  TString name = Qgr->GetName();
  TString binnedname = TString::Format("%s_binnned", name.Data()); 
  BinnedGr->SetName(binnedname);

  std::cout << "DEBUG::0.5" << std::endl;
  
  for (int i=0;i<nbins;i++)
    {

      double weighted_sum = 0;
      double weight = 0;
      
      for (int j=0;j<x.size();j++)
	{
	  if ( x[j] > binEdges[i] && x[j] < binEdges[i+1] )
	    {
	      weighted_sum += y[j]*w[j];
	      weight += w[j];
	    }
	}

      double val = weighted_sum/weight;
      BinnedGr->SetPoint(BinnedGr->GetN(), binLabels[i], val );
      
    }

  std::cout << "DEBUG::0.6" << std::endl;
  
  return BinnedGr;
  
}

double CumulantProfileContainer::Avg_M(int cent,int index)
{
  return (double) _prof_vec[index]->GetBinContent(cent);
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
  _m[1] = m_r_s[1][1]*m_r_s[1][1];
  _m[2] = m_r_s[2][1];
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

  std::map<TString,double> m;
  std::map<TString,double> q;
  
  m["m11"] = Avg_M(cent,0);
  m["m11m11"] = Avg_M(cent,1);
  m["m21"] = Avg_M(cent,2);
  m["m22"] = Avg_M(cent,3);
  m["m11m11m11"] = Avg_M(cent,4);
  m["m11m21"] = Avg_M(cent,5);
  m["m11m22"] = Avg_M(cent,6);
  m["m31"] = Avg_M(cent,7);
  m["m32"] = Avg_M(cent,8);
  m["m33"] = Avg_M(cent,9);
  m["m11m11m11m11"] = Avg_M(cent,10);
  m["m11m11m21"] = Avg_M(cent,11);
  m["m11m11m22"] = Avg_M(cent,12);
  m["m11m31"] = Avg_M(cent,13);
  m["m21m21"] = Avg_M(cent,14);
  m["m22m22"] = Avg_M(cent,15);
  m["m11m32"] = Avg_M(cent,16);
  m["m11m33"] = Avg_M(cent,17);
  m["m21m22"] = Avg_M(cent,18);
  m["m41"] = Avg_M(cent,19);
  m["m42"] = Avg_M(cent,20);
  m["m43"] = Avg_M(cent,21);
  m["m44"] = Avg_M(cent,22);

  //Convert the moments to cumulants 

  //  if ( cent > 50 && cent < 60 )
  //    {
  //      cout << "m11m11m11=" << m["m11m11m11"] << std::endl;
  //      cout << "m11m11   =" << m["m11m11"] << std::endl;
  //      cout << "m11   =" << m["m11m11"] << std::endl; 
  //    }
  
  q["q11"]         = m["m11"];
  q["q11q11"]      = m["m11m11"] - m["m11"]*m["m11"];
  q["q21"]         = m["m21"];
  q["q22"]         = m["m22"]; 
  q["q11q11q11"]   = m["m11m11m11"] - 3*m["m11m11"]*m["m11"] + 2*m["m11"]*m["m11"]*m["m11"];
  q["q11q21"]      = m["m11m21"] - 1.0*m["m11"]*m["21"];

  cout << "m11m21=" << m["m11m21"] << endl;
  cout << "  m11=" << m["m11"] << endl;
  cout << "  m21=" << m["m21"] << endl;

  cout << "q11q21=" << q["q11q21"] << endl;
  
  q["q11q22"]      = m["m11m22"] - m["m11"]*m["22"];
  q["q31"]         = m["m31"];
  q["q32"]         = m["m32"];
  q["q33"]         = m["m33"];
  q["q11q11q11q11"]= m["m11m11m11m11"] - 4*m["m11m11m11"]*m["m11"] -3*m["m11m11"]*m["m11m11"] +
                     12*m["m11m11"]*m["m11"]*m["m11"] -6*m["m11"]*m["m11"]*m["m11"]*m["m11"];
  q["q11q11q21"]   = m["m11m11m21"] -  q["q11q11"]*m["m21"] - 2*q["q11q21"]*m["m11"] - m["m11"]*m["m11"]*m["m21"];
  q["q11q11q22"]   = m["m11m11m22"] -  q["q11q11"]*m["m22"] - 2*q["q11q22"]*m["m11"] - m["m11"]*m["m11"]*m["m22"];
  q["q11q31"]      = m["m11m31"] - m["m11"]*m["m31"];
  q["q21q21"]      = m["m21m21"] - m["m21"]*m["m21"];
  q["q22q22"]      = m["m22m22"] - m["m22"]*m["m22"];
  q["q11q32"]      = m["m11m32"] - m["m11"]*m["m32"];
  q["q11q33"]      = m["m11m33"] - m["m11"]*m["m33"];
  q["q21q22"]      = m["m21m22"] - m["m21"]*m["m22"];
  q["q41"]         = m["m41"];
  q["q42"]         = m["m42"];
  q["q43"]         = m["m43"];
  q["q44"]         = m["m44"];

  return q;
   
}
