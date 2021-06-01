#include "TFile.h"
#include "TString.h"
#include "TGraphErrors.h"
#include <iostream>
#include <TAttMarker.h>
#include "TLine.h"
#include "TH1F.h"
#include "TColor.h"
#include "TCanvas.h"

using namespace std;


TCanvas *  MakeRatioCanvas();

void DrawPoint(TGraphErrors * gr, double x, double y);
void SetMarkerVals(TMarker * mark, TGraphErrors * gr);
void LoadGraphs(TString end,TFile * datafile, TGraphErrors ** gr_data_stat, TGraphErrors ** gr_data_sys);
void SetGraphVals(TGraphErrors * gr, int lineColor,double lineWidth, int lineStyle, int markerColor, double markerSize,int markerStyle,int fillColor=-999);
TGraphErrors * TranslateGraph(TGraphErrors * gr, float delta);

void DrawCumulant()
{

 gStyle->SetNumberContours(255);
 gStyle->SetOptStat(0);

  TGraphErrors * gr_pt_stat[7];
  TGraphErrors * gr_pt_base[7];
  TGraphErrors * gr_pt_sys[7];

  TGraphErrors * gr_rap_stat[7];
  TGraphErrors * gr_rap_base[7];
  TGraphErrors * gr_rap_sys[7];

  TFile * datafile = new TFile("../rootfiles/cumulants/scan_data.root","read");

  LoadGraphs( "rap", datafile, gr_rap_stat, gr_rap_sys );
  LoadGraphs(  "pt", datafile, gr_pt_stat, gr_pt_sys );

  for ( int i=0;i<7;i++)
    {
      gr_pt_base[i] = (TGraphErrors*) gr_pt_stat[i]->Clone();
      gr_rap_base[i] = (TGraphErrors*) gr_rap_stat[i]->Clone();
    }

  int turq = TColor::GetColor( 30, 132, 127);

  // C2/C1
  SetGraphVals(gr_rap_stat[4], kBlack, 3, 0, kBlack, 2.1, 21);  
  SetGraphVals(gr_pt_stat[4],  kBlack, 3, 0, kBlack, 2.1, 21);  

  SetGraphVals(gr_rap_sys[4], kGray+2, 6, 0, kBlack, 2.1, 21);  
  SetGraphVals(gr_pt_sys[4],  kGray+2, 6, 0, kBlack, 2.1, 21);  

  SetGraphVals(gr_rap_base[4], kBlack, 1, 7, kBlack, 2.1, 21);  
  SetGraphVals(gr_pt_base[4],  kBlack, 1, 7, kBlack, 2.1, 21);  


  // C3/C2
  SetGraphVals(gr_rap_stat[5], kBlack, 3, 0, kBlack, 2.1, 21);  
  SetGraphVals(gr_pt_stat[5],  kBlack, 3, 0, kBlack, 2.1, 21);  

  SetGraphVals(gr_rap_sys[5], kGray+2, 6, 0, kBlack, 2.1, 21);  
  SetGraphVals(gr_pt_sys[5],  kGray+2, 6, 0, kBlack, 2.1, 21);  

  SetGraphVals(gr_rap_base[5], kBlack, 1, 7, kBlack, 2.1, 21);  
  SetGraphVals(gr_pt_base[5],  kBlack, 1, 7, kBlack, 2.1, 21);  


  // C4/C2
  SetGraphVals(gr_rap_stat[6], kBlack, 3, 0, kBlack, 2.1, 21);  
  SetGraphVals(gr_pt_stat[6],  kBlack, 3, 0, kBlack, 2.1, 21);  

  SetGraphVals(gr_rap_sys[6], kGray+2, 6, 0, kBlack, 2.1, 21);  
  SetGraphVals(gr_pt_sys[6],  kGray+2, 6, 0, kBlack, 2.1, 21);  

  SetGraphVals(gr_rap_base[6], kBlack, 1, 7, kBlack, 2.1, 21);  
  SetGraphVals(gr_pt_base[6],  kBlack, 1, 7, kBlack, 2.1, 21);  


  TCanvas * can = new TCanvas("can","",600,900);
  TPad * pads[6];

  double x1[6];
  double x2[6];
  double y1[6];
  double y2[6];

  double ml[6];
  double mr[6];
  double mt[6];
  double mb[6];

  /******
    4 5
    2 3
    0 1
   *****/
  double mid=0.57;

  x1[0] = 0.0;
  x2[0] = mid;
  y1[0] = 0.0;
  y2[0] = 0.415;

  x1[1] = mid;
  x2[1] = 1.0;
  y1[1] = 0.0;
  y2[1] = 0.415;

  x1[2] = 0.0;
  x2[2] = mid;
  y1[2] = 0.415;
  y2[2] = 0.7;

  x1[3] = mid;
  x2[3] = 1.0;
  y1[3] = 0.415;
  y2[3] = 0.7;

  x1[4] = 0.0;
  x2[4] = mid;
  y1[4] = 0.7;
  y2[4] = 1.0;

  x1[5] = mid;
  x2[5] = 1.0;
  y1[5] = 0.7;
  y2[5] = 1.0;

  ml[0] = 0.27;
  mr[0] = 0.0;
  mt[0] = 0.0;
  mb[0] = 0.305;  

  ml[1] = 0.0;
  mr[1] = 0.35;
  mt[1] = 0.0;
  mb[1] = 0.305;  

  ml[2] = 0.27;
  mr[2] = 0.0;
  mt[2] = 0.0;
  mb[2] = 0.0;  

  ml[3] = 0.0;
  mr[3] = 0.35;
  mt[3] = 0.0;
  mb[3] = 0.0;  

  ml[4] = 0.27;
  mr[4] = 0.0;
  mt[4] = 0.05;
  mb[4] = 0.0;  

  ml[5] = 0.0;
  mr[5] = 0.35;
  mt[5] = 0.05;
  mb[5] = 0.0;  

  for (int i=0;i<6;i++)
    {
      can->cd();
      pads[i] = new TPad(TString::Format("pad_%i",i),"",x1[i],y1[i],x2[i],y2[i]);
      pads[i]->SetTopMargin(mt[i]);
      pads[i]->SetBottomMargin(mb[i]);
      pads[i]->SetLeftMargin(ml[i]);
      pads[i]->SetRightMargin(mr[i]);
      pads[i]->Draw();
    }  

  pads[0]->cd();
  gPad->SetTicks(1,1);
  //  TH1F * frame1 = gPad->DrawFrame(-0.54,-2.9,-0.14,2.9);
  TH2F * frame1 = new TH2F("2dframe","",1000,-0.54,-0.14,1000,-2.9,2.9);
  frame1->Draw();
  
  
//  frame1->GetYaxis()->SetBinLabel( frame1->GetYaxis()->FindBin(-2.0) , "-2.0");
  //  frame1->GetYaxis()->SetBinLabel( frame1->GetYaxis()->FindBin(-2.0) , "-2.0");
  
  //  TPad * coverup = new TPad("coverup", "",0.095,0.2837466,0.226,0.9421488);
  //  TPad * coverup = new TPad("coverup", "",0.135246,0.2865014,0.2677761,0.9449036);
  TPad * coverup1 = new TPad("coverup1", "",0.1391864,0.15,0.2715394-0.015,0.9421488);
  TPad * coverup = new TPad("coverup", "",0.1391864-0.005,0.2837466,0.2715394-0.005,0.9421488);
  coverup->SetFillStyle(0);
  coverup1->Draw();
  pads[0]->cd();
  coverup->Draw();
  coverup->cd();

  TLatex * newAxis = new TLatex();
  newAxis->SetTextFont(43);
  newAxis->SetTextSize(23);
  newAxis->DrawLatexNDC(0.0,0.16,"#font[122]{-}2.0");
  newAxis->DrawLatexNDC(0.275,0.527,"0.0");
  newAxis->DrawLatexNDC(0.275,0.89,"2.0");
  
  pads[0]->cd();
  //  frame1->GetYaxis()->SetTitle("Cumulant Ratios");
  frame1->GetYaxis()->SetTitle("");
  frame1->GetYaxis()->SetTitleOffset(0.65);
  frame1->GetXaxis()->SetTitleOffset(0.95);
  frame1->GetYaxis()->SetTitleSize(0.06);
  frame1->GetYaxis()->SetLabelFont(43);
  frame1->GetYaxis()->SetLabelSize(23);
  frame1->GetXaxis()->SetLabelFont(43);
  frame1->GetXaxis()->SetLabelSize(23);
  frame1->GetYaxis()->SetDecimals();

  frame1->GetYaxis()->CenterTitle();
  frame1->GetXaxis()->CenterTitle();
  frame1->GetXaxis()->SetTitle("rapidity_{min}");
  frame1->GetXaxis()->SetTitleSize(0.08);
  //  frame1->GetXaxis()->SetLabelSize(0.07);
  frame1->GetXaxis()->SetNdivisions(506);
  frame1->GetYaxis()->SetNdivisions(504);
  frame1->GetYaxis()->SetLabelOffset(0.012);

  frame1->GetXaxis()->SetTitleFont(43);
  frame1->GetXaxis()->SetTitleSize(30);
  frame1->GetXaxis()->SetTitleOffset(2.5);

  frame1->GetYaxis()->SetTitleFont(43);
  frame1->GetYaxis()->SetTitleSize(30);
  frame1->GetYaxis()->SetTitleOffset(3.4);
  frame1->GetYaxis()->SetTitle("C_{4}/C_{2}");

  pads[1]->cd();
  gPad->SetRightMargin(0.03);
  gPad->SetTicks(1,1);

  TH1F * frame2 = gPad->DrawFrame(0.6,-2.9,2.23,2.9);
  frame2->GetXaxis()->SetTitle("p^{max}_{T} (GeV/c)");
  frame2->GetXaxis()->SetTitleOffset(0.95);
  frame2->GetXaxis()->SetTitleSize(0.08);
  //  frame2->GetXaxis()->SetLabelSize(0.07);
  frame2->GetXaxis()->SetNdivisions(506);
  frame2->GetYaxis()->SetNdivisions(504);
  frame2->GetXaxis()->SetLabelFont(43);
  frame2->GetXaxis()->SetLabelSize(23);
  frame2->GetXaxis()->SetDecimals();

  frame2->GetXaxis()->SetTitleFont(43);
  frame2->GetXaxis()->SetTitleSize(30);
  frame2->GetXaxis()->SetTitleOffset(2.5);
  frame2->GetXaxis()->CenterTitle();

  pads[2]->cd();
  gPad->SetTicks(1,1);
  TH1F * frame3 = gPad->DrawFrame(-0.54,0.55,-0.14,1.35);

  //  frame3->GetYaxis()->SetTitle("Cumulant Ratios");
  frame3->GetYaxis()->SetTitle("");
  //frame3->GetYaxis()->SetTitle("");
  frame3->GetYaxis()->SetLabelFont(43);
  frame3->GetYaxis()->SetLabelSize(23);

  frame3->GetYaxis()->SetTitleFont(43);
  frame3->GetYaxis()->SetTitleSize(35);

  frame3->GetYaxis()->SetTitleOffset(0.7);
  //  frame3->GetYaxis()->SetTitleSize(0.13);
  //  frame3->GetYaxis()->SetLabelSize(0.09);
  frame3->GetYaxis()->CenterTitle();
  frame3->GetXaxis()->SetTitle("y_{min}");
  frame3->GetXaxis()->SetTitleSize(0.05);
  frame3->GetXaxis()->SetLabelSize(0.07);
  frame3->GetXaxis()->SetNdivisions(506);
  frame3->GetYaxis()->SetNdivisions(504);
  frame3->GetYaxis()->SetLabelOffset(0.012);

  frame3->GetYaxis()->SetTitleFont(43);
  frame3->GetYaxis()->SetTitleSize(30);
  frame3->GetYaxis()->SetTitleOffset(3.4);
  frame3->GetYaxis()->SetTitle("C_{3}/C_{2}");
  frame3->GetYaxis()->SetDecimals();

  pads[3]->cd();
  gPad->SetRightMargin(0.03);
  gPad->SetTicks(1,1);

  TH1F * frame4 = gPad->DrawFrame(0.6,0.55,2.23,1.35);
  frame4->GetXaxis()->SetTitle("p^{max}_{T}");
  frame4->GetXaxis()->SetTitleSize(0.07);
  frame4->GetXaxis()->SetLabelSize(0.05);
  frame4->GetXaxis()->SetNdivisions(506);
  frame4->GetYaxis()->SetNdivisions(504);

  pads[4]->cd();
  gPad->SetTicks(1,1);
  TH1F * frame5 = gPad->DrawFrame(-0.54,0.65,-0.14,1.45);

  //  frame5->GetYaxis()->SetTitle("Cumulant Ratios");
  frame5->GetYaxis()->SetLabelFont(43);
  frame5->GetYaxis()->SetLabelSize(23);

  frame5->GetYaxis()->SetTitle("");
  frame5->GetYaxis()->SetTitleOffset(0.65);
  frame5->GetYaxis()->SetTitleSize(0.06);
  //  frame5->GetYaxis()->SetLabelSize(0.08);
  frame5->GetYaxis()->CenterTitle();
  frame5->GetXaxis()->SetTitle("y_{min}");
  frame5->GetXaxis()->SetTitleSize(0.05);
  frame5->GetXaxis()->SetLabelSize(0.07);
  frame5->GetYaxis()->SetLabelOffset(0.012);
  frame5->GetXaxis()->SetNdivisions(506);
  frame5->GetYaxis()->SetNdivisions(504);

  frame5->GetYaxis()->SetTitleFont(43);
  frame5->GetYaxis()->SetTitleSize(30);
  frame5->GetYaxis()->SetTitleOffset(3.4);
  frame5->GetYaxis()->SetTitle("C_{2}/C_{1}");
  frame5->GetYaxis()->SetDecimals();

  pads[5]->cd();
  gPad->SetRightMargin(0.03);
  gPad->SetTicks(1,1);

  TH1F * frame6 = gPad->DrawFrame(0.6,0.65,2.23,1.45);
  frame6->GetXaxis()->SetTitle("p^{max}_{T}");
  frame6->GetXaxis()->SetTitleSize(0.05);
  frame6->GetXaxis()->SetLabelSize(0.07);
  frame6->GetXaxis()->SetNdivisions(506);
  frame6->GetYaxis()->SetNdivisions(504);

  //  gr_rap_base[4]->Draw("l");
  //  gr_rap_base[5]->Draw("l");

  pads[4]->cd();

  pads[2]->cd();

  pads[0]->cd();

  TLine * l1 = new TLine(-0.53,1,-0.11,1);
  l1->SetLineStyle(7);
  l1->Draw();

  TLine * l1_0 = new TLine(-0.53,0,-0.11,0);
  l1_0->SetLineStyle(2);
  l1_0->Draw();

  //  gr_rap_base[6]->Draw("l");

  //  for(int i=4;i<7;i++)
  //    {
  //      gr_rap_sys[i]->Draw("pz");
  //      gr_rap_stat[i]->Draw("pz");
  //    }

  pads[4]->cd();
  l1->Draw();
  gr_rap_base[4]->Draw("l");
  gr_rap_sys[4]->Draw("pz");
  gr_rap_stat[4]->Draw("pz");

  pads[2]->cd();
  l1->Draw();
  gr_rap_base[5]->Draw("l");
  gr_rap_sys[5]->Draw("pz");
  gr_rap_stat[5]->Draw("pz");
    
  pads[0]->cd();
  l1->Draw();
  gr_rap_base[6]->Draw("l");
  gr_rap_sys[6]->Draw("pz");
  gr_rap_stat[6]->Draw("pz");


  TLatex * txt = new TLatex();
  txt->SetTextFont(43);
  txt->SetTextSize(25);
  pads[4]->cd();
  txt->DrawLatexNDC(0.8,0.80,"(1a)");
  pads[2]->cd();

  txt->DrawLatexNDC(0.8,0.84,"(2a)");
  pads[0]->cd();

  txt->DrawLatexNDC(0.8,0.89,"(3a)");



  pads[5]->cd();
  txt->DrawLatexNDC(0.7,0.80,"(1b)");
  pads[3]->cd();

  txt->DrawLatexNDC(0.7,0.84,"(2b)");
  pads[1]->cd();

  txt->DrawLatexNDC(0.7,0.89,"(3b)");


  txt->SetTextSize(0.07);
  /*  pads[4]->cd();
  txt->DrawLatexNDC(0.65,0.75,"0.4<p_{T}<2.0 (GeV/c)");
  pads[2]->cd();
  txt->SetTextSize(0.077);
  txt->DrawLatexNDC(0.65,0.87,"0.4<p_{T}<2.0 (GeV/c)");
  */
  pads[0]->cd();
  txt->SetTextSize(0.06);
  txt->SetTextFont(43);
  txt->SetTextSize(22);
  txt->DrawLatexNDC(0.35,0.38,"0.4 < p_{T} < 2.0 (GeV/c)");
  txt->DrawLatexNDC(0.35,0.46,"y_{min} < y < 0");

  /*
  txt->SetTextSize(0.05);
  pads[5]->cd();
  txt->DrawLatexNDC(0.65,0.70,"-0.5<y<0");
  pads[3]->cd();
  txt->SetTextSize(0.05);
  txt->DrawLatexNDC(0.65,0.70,"-0.5<y<0");
  */
  pads[1]->cd();
  //  txt->SetTextSize(0.07);
  txt->DrawLatexNDC(0.11,0.38,"0.4 < p_{T} < p_{T}^{max} (GeV/c)");
  txt->DrawLatexNDC(0.11,0.46,"-0.5 < y < 0");

  txt->SetTextSize(0.04);
  txt->DrawLatex(-0.29,-0.28,"3.0 GeV");
  //  txt->DrawLatex(-0.23,-0.28,"UrQMD");

  DrawPoint( gr_rap_stat[4],-0.265,-0.47);
  DrawPoint( gr_rap_sys[4] ,-0.265,-0.47);
  DrawPoint( gr_rap_stat[5],-0.265,-0.77);
  DrawPoint( gr_rap_sys[5] ,-0.265,-0.77);
  DrawPoint( gr_rap_stat[6],-0.265,-1.06);
  DrawPoint( gr_rap_sys[6] ,-0.265,-1.06);

  TLine * lc2 = new TLine(-0.22, -0.47, -0.19, -0.47);
  TLine * lc3 = new TLine(-0.22, -0.77, -0.19, -0.77);
  TLine * lc4 = new TLine(-0.22, -1.06, -0.19, -1.06);


  lc2->SetLineWidth(12);
  lc3->SetLineWidth(12);
  lc4->SetLineWidth(12);

  lc2->Draw();
  lc3->Draw();
  lc4->Draw();

  //  can->cd(2);
  //  pads[1]->cd();

  //  gr_pt_base[4]->Draw("l");
  //  gr_pt_base[5]->Draw("l");
  pads[5]->cd();
  pads[3]->cd();
  pads[1]->cd();


  TLine * l2 = new TLine(0.6,1,2.23,1);
  l2->SetLineStyle(7);
  l2->Draw();

  TLine * l0_2 = new TLine(0.6,0,2.23,0);
  l0_2->SetLineStyle(2);
  l0_2->Draw();
  
  //  gr_pt_base[6]->Draw("l");
  
  //for(int i=4;i<7;i++)
  //    {
  //      gr_pt_sys[i]->Draw("pz");
  //      gr_pt_stat[i]->Draw("pz");
  //    }

  pads[5]->cd();
  l2->Draw();
  gr_pt_base[4]->Draw("l");
  gr_pt_sys[4]->Draw("pz");
  gr_pt_stat[4]->Draw("pz");

  pads[3]->cd();
  l2->Draw();
  gr_pt_base[5]->Draw("l");
  gr_pt_sys[5]->Draw("pz");
  gr_pt_stat[5]->Draw("pz");
    
  pads[1]->cd();
  l2->Draw();
  gr_pt_base[6]->Draw("l");
  gr_pt_sys[6]->Draw("pz");
  gr_pt_stat[6]->Draw("pz");

  //  pads[4]->cd();

  TPad * titlePad = new TPad("titlePad","",0.17,0.68,0.92,0.80);
  can->cd();
  titlePad->Draw();
  titlePad->cd();
  titlePad->SetFillStyle(0);

  //titlePad->SetFrameFillColor(0);
  //  titlePad->SetFrameFillStyle(0);
  //  titlePad->SetFrameLineColor(0);
  //  titlePad->SetFrameBorderMode(0); 

  TH1F * titleframe = (TH1F*) titlePad->DrawFrame(0,0,1,1);
  titleframe->GetXaxis()->SetNdivisions(0);
  titleframe->GetYaxis()->SetNdivisions(0);

  //  TLegend * leg = new TLegend(0.1,0.1,0.9,0.9); //new TLegend(0.25,0.05,0.87,0.25);
  //  leg->SetTextFont(42);
  //  leg->SetTextFont(43);
  //  leg->SetLegendTextSize(15);
  //  leg->SetHeader("");
  //  leg->AddEntry( gr_pt_stat[6],"Au+Au #sqrt{s_{NN}} = 3.0 GeV","p");
  //  leg->AddEntry( gr_pt_sim[6],"UrQMD","l");
  //  leg->SetLineColor(0);
  //  leg->Draw();

  TLatex * centLabel = new TLatex();
  centLabel->SetTextFont(43);
  centLabel->SetTextSize(24);

  centLabel->DrawLatexNDC(0.17,0.63,"0-5% Central Au+Au Collisions");
  centLabel->SetTextSize(20);
  centLabel->DrawLatexNDC(0.25,0.3,"#sqrt{s_{NN}} = 3.0 GeV");
  //  centLabel->DrawLatexNDC(0.67,0.3,"UrQMD");

  TGraphErrors * dataLeg = (TGraphErrors*) gr_pt_stat[6]->Clone();
  while (dataLeg->GetN()!=0) { dataLeg->RemovePoint(dataLeg->GetN()-1); }
  
  dataLeg->SetPoint(0,0.13,0.33);
  dataLeg->Draw("p");

  pads[5]->cd();

  //  TLatex * centLabel = new TLatex();
  //  centLabel->SetTextSize(0.09);
  //  centLabel->DrawLatexNDC(0.35,0.8,"Top 5% Central.");

  //  TPad * ytitle = new TPad("yTitle","",0.0,0.3,0.1,0.6);
  //  can->cd();
  //  ytitle->Draw();
  //  TLatex * title_latex = new TLatex();
  //  title_latex->SetTextAngle(90);
  //  title_latex->DrawLatexNDC(0.5,0.5,"Cumulant Ratios");

  
  can->SaveAs("../img/CumulantRatios.png");

}

void DrawPoint(TGraphErrors * gr, double x, double y)
{

  TGraphErrors * gr_point = (TGraphErrors*) gr->Clone();

  while ( gr_point->GetN() != 0 )
    {
      gr_point->RemovePoint(gr_point->GetN()-1);
    }

  gr_point->SetMarkerSize(3.2);
  gr_point->SetPoint(0,x,y);
  gr_point->Draw("p");

}


void LoadGraphs(TString end,TFile * datafile, TGraphErrors ** gr_data_stat, TGraphErrors ** gr_data_sys)
{
  TString gr_name_stat[7] = {"c1","c2","c3","c4","cr21","cr32","cr42"};

  for (int i=0;i<7;i++)
    {
      gr_data_stat[i] = (TGraphErrors*) datafile->Get( TString::Format("%s_%s_stat",end.Data(),gr_name_stat[i].Data() ) );
      gr_data_sys[i]  = (TGraphErrors*) datafile->Get( TString::Format("%s_%s_sys",end.Data(),gr_name_stat[i].Data() ) );
    }

}


void SetGraphVals(TGraphErrors * gr, 
		  int lineColor,   double lineWidth, int lineStyle, 
		  int markerColor, double markerSize,int markerStyle, int fillColor)
{

  for (int i=0;i<7;i++)
    {
      gr->SetLineColor(lineColor);
      gr->SetLineWidth(lineWidth);
      gr->SetLineStyle(lineStyle);

      gr->SetMarkerColor(markerColor);
      gr->SetMarkerStyle(markerStyle);
      gr->SetMarkerSize(markerSize);

       if ( fillColor != -999 ) gr->SetFillColor( fillColor );
    }
}

 
TGraphErrors * TranslateGraph(TGraphErrors * gr, float delta)
{
  TGraphErrors * t_gr = new TGraphErrors();
  for ( int i=0;i<gr->GetN();i++)
    {
      t_gr->SetPoint(t_gr->GetN(),gr->GetX()[i]+delta,gr->GetY()[i]);
      t_gr->SetPointError(t_gr->GetN()-1,gr->GetEX()[i],gr->GetEY()[i]);
    }
  return t_gr;
}
