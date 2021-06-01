

void GetTop5(TGraphErrors * gr,double &y,double &yerr)
{
  for (int i=0;i<gr->GetN();i++)
    {
      if ( gr->GetX()[i] > 320 )
	{
	  y = gr->GetY()[i];
	  yerr = gr->GetEY()[i];
	  return;
	}
    }
}

void MakeTopFive()
{

  TFile * f_rp[4];
  TFile * f_pt[4];

  TGraphErrors * g_rp[2][7][4];
  TGraphErrors * g_pt[2][7][4];

  TGraphErrors * g5_rp[2][7];
  TGraphErrors * g5_pt[2][7];

  TString gr_name_stat[7] = {"c1","c2","c3","c4","cr21","cr32","cr42"};
  TString gr_name_sys[7] = {"c1_sys","c2_sys","c3_sys","c4_sys","cr21_sys","cr32_sys","cr42_sys"};

  //Load files
  for (int i=0;i<4;i++)
    {
      f_rp[i] = new TFile(TString::Format("../rootfiles/cumulants/cumulants_n0p%i_0.root",i+2 ),"read");
      if ( i <3 ) f_pt[i] = new TFile(TString::Format("../rootfiles/cumulants/cumulants_n0p5_0_pt%i.root",i+1 ),"read");
    }

  f_pt[3] = f_rp[3];

  TCanvas * can = new TCanvas("can","",1100,800);
  can->Divide(6,4);
  int canIndex=1;
  
  //Load Graphs
  for (int i=0;i<7;i++)
    {
      for (int ii=0;ii<4;ii++)
	{
	  can->cd(canIndex); canIndex++;
	  //STAT
	  g_rp[0][i][ii] = (TGraphErrors*) f_rp[ii]->Get(gr_name_stat[i]);
	  g_pt[0][i][ii] = (TGraphErrors*) f_pt[ii]->Get(gr_name_stat[i]);
	  //SYS
	  g_rp[1][i][ii] = (TGraphErrors*) f_rp[ii]->Get(gr_name_sys[i]);
	  g_pt[1][i][ii] = (TGraphErrors*) f_pt[ii]->Get(gr_name_sys[i]);
	  
	  g_rp[1][i][ii]->SetLineColor(kRed);
	  g_pt[1][i][ii]->SetLineColor(kBlue);

	  g_rp[0][i][ii]->SetMarkerStyle(kRed);
	  g_pt[0][i][ii]->SetMarkerStyle(kBlue);

	  g_rp[1][i][ii]->Draw("ap");
	  g_pt[1][i][ii]->Draw("lp");

	  g_rp[0][i][ii]->Draw("lp");
	  g_pt[0][i][ii]->Draw("lp");
	  
	}
    }
    

  double pt_vals[4] = {0.8, 1.2, 1.6, 2.0};
  double rap_vals[4] = {-0.2, -0.3, -0.4, -0.5};

  for (int i=0;i<2;i++)
    {
      for (int ii=0;ii<7;ii++)
	{
	  g5_rp[i][ii] = new TGraphErrors();
	  g5_rp[i][ii]->SetName(TString::Format("rap_%s_%s",gr_name_stat[ii].Data(), i ? "sys":"stat"));
	  g5_pt[i][ii] = new TGraphErrors();
	  g5_pt[i][ii]->SetName(TString::Format("pt_%s_%s",gr_name_stat[ii].Data(), i ? "sys":"stat"));

	  for (int iii=0;iii<4;iii++)
	    {
	      double y_rap, y_rap_err;
	      double y_pt, y_pt_err;

	      GetTop5( g_rp[i][ii][iii], y_rap, y_rap_err );
	      GetTop5( g_pt[i][ii][iii], y_pt, y_pt_err );

	      g5_rp[i][ii]->SetPoint( g5_rp[i][ii]->GetN(), rap_vals[iii], y_rap );
	      g5_rp[i][ii]->SetPointError( g5_rp[i][ii]->GetN()-1, 0, y_rap_err );
	      
	      g5_pt[i][ii]->SetPoint( g5_pt[i][ii]->GetN(), pt_vals[iii], y_pt );
	      g5_pt[i][ii]->SetPointError( g5_pt[i][ii]->GetN()-1, 0, y_pt_err );
	    }
	}
    }


  TFile * fout = new TFile("../rootfiles/cumulants/scan_data.root","recreate");
  fout->cd();

  for (int i=0;i<2;i++)
    {
      for (int ii=0;ii<7;ii++)
	{
	  g5_rp[i][ii]->Write();
	  g5_pt[i][ii]->Write();
	}
    }

  
}
