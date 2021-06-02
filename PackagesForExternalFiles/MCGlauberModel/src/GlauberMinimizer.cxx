#include "GlauberMinimizer.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph2D.h"
#include "GlauberClass.h"
#include "GlauberUtil.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TString.h"

GlauberMinimizer::GlauberMinimizer(TH1D* a_measuredMultiplicity, TFile* a_glauberFile, double a_lowMultCut, double a_highMultCut){
  m_multiplicity = (TH1D*) a_measuredMultiplicity->Clone();
  m_multiplicity->SetName("dataHisto");
  m_simMult = (TH1D*) m_multiplicity->Clone();
  m_simMult->SetName("simHisto");
  m_glauberEvent = new GlauberClass();


  TTree* glauberTree = (TTree*) a_glauberFile->Get("GlauberTree");
  glauberTree->FindBranch("GlauberData")->SetAddress(&m_glauberEvent);

  m_numGlauberEntries = glauberTree->GetEntries();
  for(int iii = 0; iii < m_numGlauberEntries; iii++){
    glauberTree->GetEntry(iii);
    m_nPart.push_back(m_glauberEvent->GetNParticipants());
    m_nColl.push_back(m_glauberEvent->GetNBinaryCollisions());
  }//End Loop Over Glauber Tree

  m_lowMultCut  = a_lowMultCut;
  m_highMultCut = a_highMultCut;
  m_lowMultCutBin  = m_multiplicity->FindBin(m_lowMultCut);
  m_highMultCutBin = m_multiplicity->FindBin(m_highMultCut);
  m_multiplicity->Scale( 1.0/ ((double)m_multiplicity->Integral(m_lowMultCutBin,m_highMultCutBin)) );  

  m_chiSqrdGraph = new TGraph2D();
  m_chiSqrdGraph->SetName("chiSqrdGraph");
  m_glauberEvent = new GlauberClass();

  m_chiSqrdMethod = "normal";
}

void GlauberMinimizer::FitGlauberModel( double a_initialNpp, double a_initialK, double a_initialX, int a_numEventOverride, TString outdir, TString ofile){


  m_minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  ROOT::Math::Functor* chiSquaredFunctor = new ROOT::Math::Functor(this,&GlauberMinimizer::chiSquaredFunct,3);
  m_minimizer->SetFunction(*chiSquaredFunctor);
  if(a_numEventOverride > 0){
    m_numGlauberEntries = a_numEventOverride;
  }

  double stepSize = 10;

  m_minimizer->SetVariable(0, "Npp" ,a_initialNpp, stepSize);
  m_minimizer->SetVariableLimits(0,0.1, 2.0);
  //m_minimizer->SetVariable(0.0,10.0 * avgTemp);

  m_minimizer->SetVariable(1,"k", a_initialK, stepSize);
  m_minimizer->SetVariableLimits(1, 0.0, 5.0 );

  if(a_initialX < 0.0){
    m_minimizer->SetVariable(2,"x", 0.12, stepSize);
    m_minimizer->FixVariable(2);
  }else{
    m_minimizer->SetVariable(2,"x", a_initialX, stepSize);
    m_minimizer->SetVariableLimits(2, 0.0, 1.0 );
  }

  m_minimizer->SetPrintLevel(5);
  m_minimizer->SetMaxIterations(50000);
  m_minimizer->SetMaxFunctionCalls(1E9);
  m_minimizer->SetTolerance(0.001);
  //m_minimizer->SetValidError(true);
  m_minimizer->Minimize();
  int status = m_minimizer->Status();
  cout << "Fit Status: " << status << endl;
  if(status == 0){
    cout << "Good Fit" << endl;
  }
  else{
    cout << "Bad Fit" << endl;
    //continue;
  }

  TCanvas* finalCanvas = new TCanvas("canvas");
  finalCanvas->cd();
  gPad->SetRightMargin(0.4);
  m_multiplicity->SetLineColor(kBlue);
  m_simMult->SetLineColor(kRed);
  //m_multiplicity->Draw("E");
  m_multiplicity->Draw();
  m_simMult->Draw("same");
  std::stringstream namestring1;
  namestring1.str("");
  namestring1<<outdir<<"finalFit.png";
  finalCanvas->SaveAs(namestring1.str().c_str());

  //TFile* outFile = new TFile("fitResults.root","RECREATE");
  TFile* outFile = new TFile(ofile,"RECREATE");
  m_simMult->Write();
  m_multiplicity->Write();
  m_chiSqrdGraph->Write();

  for(int iii = 0; iii < (int) m_simHistory.size(); iii++){
    m_simHistory[iii]->Write();
  }



  //image for search space
  TCanvas* canvHisotry = new TCanvas("canvasHistory");
  canvHisotry->cd();
  gPad->SetRightMargin(0.4);
  gPad->SetLogy();
  gPad->SetRightMargin(0.05);
  m_multiplicity->Draw();
  for(int iii = 0; iii < (int) m_simHistory.size(); iii++){
    m_simHistory[iii]->SetLineColor(kRed);
    m_simHistory[iii]->Draw("same");
  }
  m_multiplicity->Draw("same");
  std::stringstream namestring2;
  namestring2.str("");
  namestring2<<outdir<<"fitHistory.png";
  canvHisotry->SaveAs(namestring2.str().c_str());

  //best CHi^2 image
  double lowestChiSqrd = 1e200;
  vector<double> paramLowChi;
  int lowestIndex = -1;
  for(int iii = 0; iii < (int) m_simParamHistory.size(); iii++){
    if(m_simParamHistory[iii][3] < lowestChiSqrd){
      lowestChiSqrd = m_simParamHistory[iii][3];
      paramLowChi = m_simParamHistory[iii];
      lowestIndex = iii;
    }
  }

  TCanvas* canvLowest = new TCanvas("canvasLowest");
  canvLowest->cd();
  m_simHistory[lowestIndex]->Draw("HIST");
  m_multiplicity->Draw("same");
  TPaveText* pt = new TPaveText(.85,.8,.95,.95, "NDC");
  pt->AddText(Form("N_{pp} = %1.3f", paramLowChi[0] ));
  pt->AddText(Form("k = %1.3f", paramLowChi[1] ));
  pt->AddText(Form("x = %1.3f", paramLowChi[2] ));
  pt->AddText(Form("#chi^{2}/ndf = %1.3f", paramLowChi[3]/((double)m_multiplicity->GetNbinsX() + 3) ));
  pt->Draw();
  std::stringstream namestring3;
  namestring3.str("");
  namestring3<<outdir<<"fitBestChi.png";
  canvLowest->SaveAs(namestring3.str().c_str());





}





double GlauberMinimizer::chiSquaredFunct(const double* a_param){
  // [0] is npp  [1] is k   [2] is hardness
  //int numGlauberEvents = 1000;

  m_glauberEvent->SetNegativeBinomialParameters(a_param[0],a_param[1]);
  TH1D * nbdHist = MakeNegativeBinomialHist(a_param[0],a_param[1],a_param[2],0,*m_glauberEvent);
  for(int bin = 1; bin <= m_simMult->GetNbinsX(); bin++){
    m_simMult->SetBinContent(bin,0.0);
    m_simMult->SetBinError(bin,0.0);
  }

  //Generate the multiplicity distribution
  for(int iii = 0; iii < m_numGlauberEntries; iii++){
    m_simMult->Fill(m_glauberEvent->ProduceParticles(m_nPart[iii],m_nColl[iii],nbdHist,a_param[2],false) );
  }

  delete nbdHist;

  m_simMult->Scale( 1.0/ ((double)m_simMult->Integral(m_lowMultCutBin,m_simMult->GetNbinsX())) );


  double chiSqrd = 0.0;

  if(m_chiSqrdMethod == "normal"){

    for(int bin = m_lowMultCutBin; bin <= m_multiplicity->GetNbinsX(); bin++){
      if(m_multiplicity->GetBinError(bin) > 0.0){
        chiSqrd += ( m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin))
                     *(m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin))
                     /(m_multiplicity->GetBinError(bin)*m_multiplicity->GetBinError(bin));
      }
    }

  }else if(m_chiSqrdMethod == "geometric"){
    for(int bin = m_lowMultCutBin; bin <= m_multiplicity->GetNbinsX(); bin++){
      chiSqrd += ( m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin))
                     *(m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin));
    }
  }else if(m_chiSqrdMethod == "mixed"){

    for(int bin = m_lowMultCutBin; bin <= m_multiplicity->GetNbinsX(); bin++){
      if(m_multiplicity->GetBinError(bin) > 0.0 && m_simMult->GetBinError(bin) > 0.0){
        chiSqrd += ( m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin))
                     *(m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin))
                     /(m_multiplicity->GetBinError(bin)*m_simMult->GetBinError(bin));
      }else if(m_multiplicity->GetBinError(bin) <= 0.0 && m_simMult->GetBinError(bin) > 0.0){
        chiSqrd += ( m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin))
                     *(m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin))
                     /(m_simMult->GetBinError(bin)*m_simMult->GetBinError(bin));
      }else if(m_multiplicity->GetBinError(bin) > 0.0 && m_simMult->GetBinError(bin) <= 0.0){
        chiSqrd += ( m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin))
                     *(m_multiplicity->GetBinContent(bin) - m_simMult->GetBinContent(bin))
                     /(m_multiplicity->GetBinError(bin)*m_multiplicity->GetBinError(bin));
      }
    }

  }else if(m_chiSqrdMethod == "blah..."){

  }


  cout << "Npp: " << a_param[0] << " k: " << a_param[1] << " x: " << a_param[2] << "  chiSqrd: " << chiSqrd << endl;
  m_chiSqrdGraph->SetPoint(m_chiSqrdGraph->GetN(),a_param[0],a_param[1],1.0/chiSqrd);
  
  m_simHistory.push_back((TH1D*) m_simMult->Clone());
  //vector<double> paramChi = {a_param[0], a_param[1] , a_param[2] , chiSqrd};
  vector<double> paramChi;
  paramChi.push_back(a_param[0]);
  paramChi.push_back(a_param[1]);
  paramChi.push_back(a_param[2]);
  paramChi.push_back(chiSqrd);
  m_simParamHistory.push_back(paramChi);
  return chiSqrd;

}











