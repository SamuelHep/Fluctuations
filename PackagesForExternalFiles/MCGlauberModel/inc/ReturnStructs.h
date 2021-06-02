//These are various structures used as return trypes for functions

class TH1D;
class TH3D;
class TGraph;
class TGraph2D;
class TCanvas;
class TNtuple;

//____________________________________________
struct NegBinomialSearchResults{

  //This Structure is returned from FindBestFitNegativeBinomialParameters()

  //These are variables used as part of FindBestFitNegativeBinomialParameters()
  Double_t npp;
  Double_t k;
  Double_t hardness;
  Double_t InverseChi2;
  TH1D *bestFitHisto;
  TH1D *Prob1D;
  TGraph *chi2Graph1D;
  TGraph2D *chi2Graph2D;
  TH3D *chi3D;
  //TGraph2D *chi2Graph2Dhard0;
  //TGraph2D *chi2Graph2Dhard1;
  //TGraph2D *chi2Graph2Dhard2;
  //TGraph2D *chi2Graph2Dhard3;
  //TGraph2D *chi2Graph2Dhard4;
  //TGraph2D *chi2Graph2Dhard5;
  //TGraph2D *chi2Graph2Dhard6;
  //TGraph2D *chi2Graph2Dhard7;
  //TGraph2D *chi2Graph2Dhard8;
  //TGraph2D *chi2Graph2Dhard9; 
  TNtuple *nTup3D;

  //These are variables used as part of FindCentralityBinCuts() and FindNpartNcollDistributions()
  Int_t nCentralityBins;
  std::vector<double> centralityBinDefinitions;
  std::vector<double> centralityBinCuts;
  TH1D *nPartTotalHisto;
  TH1D *nCollTotalHisto;
  TH1D *impactParamTotalHisto;
  TH1D *nPartHistos;
  TH1D *nCollHistos;
  TH1D *impactParamHistos;

  std::vector<double> nPartMeans;
  std::vector<double> nPartStatErrors;
  std::vector<double> nPartSysErrors;
  std::vector<double> nCollMeans;
  std::vector<double> nCollStatErrors;
  std::vector<double> nCollSysErrors;
  std::vector<double> impactParamMeans;
  std::vector<double> impactParamStatErrors;
  std::vector<double> impactParamSysErrors;

  //Canvas that could be useful for future drawing
  TCanvas *nPartTotalSysErrCanvas;
  TCanvas *nCollTotalSysErrCanvas;
  TCanvas *impactParamTotalSysErrCanvas;
  TCanvas *resultsTableCanvas;
};

//______________________________________________
//This class get written to a seperate "Correction Info File" that is designed to contain
//the minimal amount of things necessary to perform the correction operations during
//looping over data.
class RefMultCorrInfo {

 public:
  Long64_t triggerID;
  Double_t refMultCentBinCuts[16];

  RefMultCorrInfo();
  ~RefMultCorrInfo();

  ClassDef(RefMultCorrInfo,1);
};
