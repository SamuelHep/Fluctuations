#include <vector>
#include <utility>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TString.h"
#include <iostream>

#include "SystematicCalculation.h"

using namespace std;

ClassImp( SystematicCalculation );

SystematicCalculation::SystematicCalculation()
{

  std::vector<TString> grnames = {"C1_cbwc","C2_cbwc","C3_cbwc","C4_cbwc","C5_cbwc","C6_cbwc","C2_C1_cbwc","C3_C2_cbwc","C4_C2_cbwc","C5_C1_cbwc","C6_C2_cbwc",
				  "K1_cbwc","K2_cbwc","K3_cbwc","K4_cbwc","K5_cbwc","K6_cbwc","K2_K1_cbwc","K3_K1_cbwc","K4_K1_cbwc","K5_K1_cbwc","K6_K1_cbwc"};


  std::vector<TString> simple_names_norm = {"c1","c2","c3","c4","c5","c6","cr21","cr32","cr42","cr51","cr62",
					    "k1","k2","k3","k4","k5","k6","kr21","kr31","kr41","kr51","kr61"};
  std::vector<TString> simple_names_sys  = {"c1_sys","c2_sys","c3_sys","c4_sys","c5_sys","c6_sys","cr21_sys","cr32_sys","cr42_sys","cr51_sys","cr62_sys",
					    "k1_sys","k2_sys","k3_sys","k4_sys","k5_sys","k6_sys","kr21_sys","kr31_sys","kr41_sys","kr51_sys","kr61_sys"};
 
  _simple_names_norm = simple_names_norm;
  _simple_names_sys  = simple_names_sys;

  _grnames = grnames;

}


SystematicCalculation::~SystematicCalculation()
{

  for ( auto &g : _nominal_gr_vec )
    {
      if ( g ) delete g;
    }

  for ( auto &g : _sys_gr_vec )
    {
      if ( g ) delete g;
    }

  for ( auto &g_vec : _nsys_single_gr )
    {
      for ( auto &g : g_vec )
	{
	  if ( g ) delete g;
	}
    }

  for ( auto &g_vec : _nsys_double_gr )
    {
      for ( auto &g : g_vec )
	{
	  if ( g ) delete g;
	}
    }
  
}

bool SystematicCalculation::AddNominal(TString filename)
{

  TFile * f = new TFile(filename,"read");
  
  vector<TGraphErrors*>  gr_norm_vec(_grnames.size(),(TGraphErrors*)NULL);
  vector<TGraphErrors*>  gr_sys_vec(_grnames.size(),(TGraphErrors*)NULL);
  int i=0;

  for ( auto &gr : gr_norm_vec )
    {
      gr = (TGraphErrors*) ((TGraphErrors*) f->Get(_grnames[i]))->Clone();
      if ( !gr ) 
	{
	  cout << "The graph " << _grnames[i] << " in the file " << filename << " was not found"  << endl;
	  return false;
	}
      i++;
    }

  i=0;

  for ( auto &gr : gr_sys_vec )
    {
      gr = (TGraphErrors*) ((TGraphErrors*) f->Get(_grnames[i]))->Clone();
      if ( !gr ) 
	{
	  cout << "The graph " << _grnames[i] << " in the file " << filename << " was not found"  << endl;
	  return false;
	}
      i++;
    }


  _nominal_gr_vec = gr_norm_vec;
  _sys_gr_vec = gr_sys_vec;

  f->Close();

  return true;

}

bool SystematicCalculation::AddSysSingle(TString filename)
{

  TFile * f = new TFile(filename,"read");
  
  vector<TGraphErrors*>  gr_vec(_grnames.size(),(TGraphErrors*)NULL);
  int i=0;

  cout << "File name = " << filename << endl;

  for ( auto &gr : gr_vec )
    {
      gr = (TGraphErrors*) ((TGraphErrors*) f->Get(_grnames[i]))->Clone();
      if ( !gr ) 
	{
	  cout << "The graph " << _grnames[i] << " in the file " << filename << " was not found"  << endl;
	  return false;
	}
      i++;
    }

  _nsys_single_gr.push_back( gr_vec );

  f->Close();

  return true;

}

bool SystematicCalculation::AddSysPair(TString filename1,TString filename2)
{

  TFile * f1 = new TFile(filename1,"read");
  TFile * f2 = new TFile(filename2,"read");
  
  vector<TGraphErrors*>  gr_vec1(_grnames.size(),(TGraphErrors*)NULL);
  vector<TGraphErrors*>  gr_vec2(_grnames.size(),(TGraphErrors*)NULL);

  cout << "File name 1 = " << filename1 << endl;
  cout << "File name 2 = " << filename2 << endl;
      
  int i=0;
  for ( auto &gr : gr_vec1 )
    {

      gr  = (TGraphErrors*) ((TGraphErrors*) f1->Get(_grnames[i]))->Clone();
      if ( !gr ) 
	{
	  cout << "The graph " << _grnames[i] << " in the file " << filename1 << " was not found"  << endl;
	  return false;
	}
      i++;
    }


  i=0;
  for ( auto &gr : gr_vec2 )
    {
      gr  = (TGraphErrors*) ((TGraphErrors*) f2->Get(_grnames[i]))->Clone();
      if ( !gr ) 
	{
	  cout << "The graph " << _grnames[i] << " in the file " << filename2 << " was not found"  << endl;
	  return false;
	}
      i++;
    }


  _nsys_double_gr.push_back( gr_vec1 );
  _nsys_double_gr.push_back( gr_vec2 );

  f1->Close();
  f2->Close();
  
  return true;

}

bool SystematicCalculation::Calculate()
{

  if (!CheckGraphVectors()) return false;
  if ( _nominal_gr_vec.size() == 0 ) return false;

  for ( int iGr=0; iGr<_nominal_gr_vec.size();iGr++)
    {
      
      //Get the nominal g
      TGraphErrors * nom_gr = _nominal_gr_vec[iGr];
      TGraphErrors * sys_gr = _sys_gr_vec[iGr];
      
     //Loop through points
      for ( int iPt=0; iPt<nom_gr->GetN(); iPt++ )
	{
	  double nom_x = nom_gr->GetX()[iPt];
	  double nom_y = nom_gr->GetY()[iPt];
	  
	  double sys_err2=0;
	  
	  for ( int iSys=0;iSys<_nsys_single_gr.size();iSys++ )
	    {
	      double sys_y = _nsys_single_gr[iSys][iGr]->GetY()[iPt];
	      cout << "single sys_y=" << sys_y << " mag=" << sqrt( pow( nom_y - sys_y, 2 )) << endl;
	      sys_err2+= pow( nom_y - sys_y, 2 );
	    }

	  for ( int iSys=0;iSys<_nsys_double_gr.size();iSys++ )
	    {
	      double sys_y_double = _nsys_double_gr[iSys][iGr]->GetY()[iPt];
	      cout << "double sys_y=" << " mag=" << sqrt(pow( nom_y - sys_y_double, 2 )) << endl;
	      sys_err2+= 0.5*(pow( nom_y - sys_y_double, 2 ));
	    }

	  double sys_err = sqrt( sys_err2 );
	  cout << "nom x=" << nom_x << " nom_y=" << nom_y <<" sys error=" << sys_err << endl;

	  sys_gr->SetPointError(iPt,0,sys_err);
	}
    }
  
  return true;

}

bool SystematicCalculation::CheckGraphVectors()
{
  int Ngr= _nominal_gr_vec.size();
  int N  = _nominal_gr_vec[0]->GetN();

  for ( auto &gr : _nominal_gr_vec )
    {
      int n = gr->GetN();
      if ( N==n ) N=n;
      else { 
	cout << "SystematicCalculation::Calculate()"
	     << "Graph->GetN() have unequal length... exit calculation" << endl; 
	return false;
      }
    }  

  for ( auto &gr_vec : _nsys_single_gr )
    {
      int ngr = gr_vec.size();
      if ( Ngr==-1 || Ngr==ngr ) Ngr=ngr;
      else { 
	cout << "SystematicCalculation::Calculate()"
	     << "GraphVectors have unequal length... exit calculation" << endl; 
	break;
      }

      for ( auto &gr : gr_vec )
	{
	  int n = gr->GetN();
	  if ( N==-1 || N==n ) N=n;
	  else { 
	    cout << "SystematicCalculation::Calculate()"
		 << "Graph->GetN() have unequal length... exit calculation" << endl; 
	    return false;
	  }
	}
    }

  for ( auto &gr_vec : _nsys_double_gr )
    {

      int ngr = gr_vec.size();
      if ( Ngr==-1 || Ngr==ngr ) Ngr=ngr;
      else { 
	cout << "SystematicCalculation::Calculate()"
	     << "GraphVectors have unequal length... exit calculation" << endl; 
	return false;
      }

      for ( auto &gr : gr_vec )
	{
	  int n = gr->GetN();
	  
	  if ( N==-1 || N==n ) N=n;
	  else { 
	    cout << "SystematicCalculation::Calculate()"
		 << "Graph->GetN() have unequal length... exit calculation" << endl; 
	    break;
	  }
	  
	}
    }

  return true;

}


bool SystematicCalculation::WriteToOutFile(TString outfilename, bool rename)
{
  
  TFile * f= new TFile( outfilename, "recreate");
  
  for ( int igr=0; igr<_nominal_gr_vec.size(); igr++)
    {

      if ( rename == true )
	{
	  _nominal_gr_vec[igr]->SetName( _simple_names_norm[igr] );
	  _sys_gr_vec[igr]->SetName( _simple_names_sys[igr] );
	}
      
      f->cd();
      _nominal_gr_vec[igr]->Write();
      _sys_gr_vec[igr]->Write();
    }
  
  f->Close();

}

