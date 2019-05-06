#ifndef ElectronEfficiencyAnalyzer_H
#define ElectronEfficiencyAnalyzer_H
 
#include "TH1.h"
#include "TH2.h"
#include "PDAnalysis/Ntu/interface/MGBaseAnalyzer.h"
#include "PDAnalysis/Ntu/interface/MGSelector.h"
#include "PDAnalysis/Ntu/interface/MGGenTools.h"
#include "MGTools/PlottingTools/interface/PlottingTools.h"

class ElectronEfficiencyAnalyzer: 
    public MGBaseAnalyzer,
    public virtual MGSelector,
    public virtual MGGenTools, 
    public virtual PlottingTools 
{
 public:

  ElectronEfficiencyAnalyzer() : MGBaseAnalyzer("ElectronEfficiencyAnalyzer") {};
  virtual ~ElectronEfficiencyAnalyzer();

  // function called before starting the analysis
  virtual void beginJob();

  // functions to book the histograms
  void book();

  // functions called for each event
  // function to reset class content before reading from file
  virtual void reset();
  // function to do event-by-event analysis,
  // return value "true" for accepted events
  virtual bool analyze( int entry, int event_file, int event_tot );

  // function called at the end of the analysis
  virtual void endJob();

  bool verbose;

 private:
  std::map<int, TH1D*> hNGenEle;
  std::map<int, TCanvas*> cNGenEle;
  std::map<int, TH1D*> hGenEleCharge;
  std::map<int, TCanvas*> cGenEleCharge;
  std::map<int, TH1D*> hGenElePt;
  std::map<int, TCanvas*> cGenElePt;
  std::map<int, TH1D*> hGenElePhi;
  std::map<int, TCanvas*> cGenElePhi;
  std::map<int, TH1D*> hGenEleEta;
  std::map<int, TCanvas*> cGenEleEta;
  std::map<int, TH2D*> hGenElePtVsEta;
  std::map<int, TCanvas*> cGenElePtVsEta;
  std::map<int, TH1D*> hGenEleE;
  std::map<int, TCanvas*> cGenEleE;
  std::map<int, TH1D*> hGenEleDb;
  std::map<int, TCanvas*> cGenEleDb;
  std::map<int, TH1D*> hGenEleBsChargeCorr;
  std::map<int, TCanvas*> cGenEleBsChargeCorr;
  std::map<int, TH1D*> hGenEleRecoChargeCorr;
  std::map<int, TCanvas*> cGenEleRecoChargeCorr;
  
  std::map<int, TH1D*> hNEle;
  std::map<int, TCanvas*> cNEle;
  std::map<int, TH1D*> hEleCharge;
  std::map<int, TCanvas*> cEleCharge;
  std::map<int, TH1D*> hElePt;
  std::map<int, TCanvas*> cElePt;
  std::map<int, TH1D*> hElePhi;
  std::map<int, TCanvas*> cElePhi;
  std::map<int, TH1D*> hEleEta;
  std::map<int, TCanvas*> cEleEta;
  std::map<int, TH2D*> hElePtVsEta;
  std::map<int, TCanvas*> cElePtVsEta;  
  std::map<int, TH1D*> hEleE;
  std::map<int, TCanvas*> cEleE;
  std::map<int, TH1D*> hEleIDs;
  std::map<int, TCanvas*> cEleIDs;
  std::map<int, TH1D*> hEleDb;
  std::map<int, TCanvas*> cEleDb;
  std::map<int, TH1D*> hEleBsChargeCorr;
  std::map<int, TCanvas*> cEleBsChargeCorr;
  std::map<int, TH1D*> hEleGenChargeCorr;
  std::map<int, TCanvas*> cEleGenChargeCorr;
  std::map<int, TH1D*> hElePassesHZZV1ID;
  std::map<int, TCanvas*> cElePassesHZZV1ID;
  std::map<int, TH1D*> hElePassesHZZV2ID;
  std::map<int, TCanvas*> cElePassesHZZV2ID;
  std::map<int, TH1D*> hElePassesMVANILID;
  std::map<int, TCanvas*> cElePassesMVANILID;
  
  std::map<int, TH1D*> hElePVZDistanceFromBestPV;
  std::map<int, TCanvas*> cElePVDistanceFromBestPV;
  
  
  std::map<int, std::vector<TH1D*>* > vhEleEtaVsPt;
  std::map<int, std::vector<TCanvas*>* > vcEleEtaVsPt;
  
  std::map<int, std::vector<TH1D*>*> vhEleHZZMVAOutputVsPt;
  std::map<int, std::vector<TCanvas*>*> vcEleHZZMVAOutputVsPt;
  std::map<int, std::vector<TH1D*>*> vhEleHZZMVAOutputVsHZZMVACat;
  std::map<int, std::vector<TCanvas*>*> vcEleHZZMVAOutputVsHZZMVACat;
  std::map<int, std::vector<TH1D*>*> vhEleHZZMVACatVsPt;
  std::map<int, std::vector<TCanvas*>*> vcEleHZZMVACatVsPt;
  
  std::map<int, std::vector<TH1D*>*> vhEleNIMVAOutputVsPt;
  std::map<int, std::vector<TCanvas*>*> vcEleNIMVAOutputVsPt;
  std::map<int, std::vector<TH1D*>*> vhEleNIMVAOutputVsNIMVACat;
  std::map<int, std::vector<TCanvas*>*> vcEleNIMVAOutputVsNIMVACat;
  std::map<int, std::vector<TH1D*>*> vhEleNIMVACatVsPt;
  std::map<int, std::vector<TCanvas*>*> vcEleNIMVACatVsPt;
  
  std::map<int, std::vector<TH1D*>*> vhEleIMVAOutputVsPt;
  std::map<int, std::vector<TCanvas*>*> vcEleIMVAOutputVsPt;
  std::map<int, std::vector<TH1D*>*> vhEleIMVAOutputVsIMVACat;
  std::map<int, std::vector<TCanvas*>*> vcEleIMVAOutputVsIMVACat;
  std::map<int, std::vector<TH1D*>*> vhEleIMVACatVsPt;
  std::map<int, std::vector<TCanvas*>*> vcEleIMVACatVsPt;

  std::map<int, TGraphAsymmErrors*> gGenEleRecoEffVsCharge;
  std::map<int, TGraphAsymmErrors*> gGenEleRecoEffVsPt;
  std::map<int, TGraphAsymmErrors*> gGenEleRecoEffVsPhi;
  std::map<int, TGraphAsymmErrors*> gGenEleRecoEffVsEta;
  std::map<int, TGraphAsymmErrors*> gGenEleRecoEffVsE;
  std::map<int, TGraphAsymmErrors*> gGenEleRecoEffVsDb;
  std::map<int, TGraphAsymmErrors*> gGenEleRecoEffVsBsChargeCorr;
  std::map<int, TGraphAsymmErrors*> gGenEleRecoEffVsRecoChargeCorr;

  std::map<int, TGraphAsymmErrors*> gEleEffVsCharge;
  std::map<int, TGraphAsymmErrors*> gEleEffVsPt;
  std::map<int, TGraphAsymmErrors*> gEleEffVsPhi;
  std::map<int, TGraphAsymmErrors*> gEleEffVsEta;
  std::map<int, TGraphAsymmErrors*> gEleEffVsE;
  std::map<int, TGraphAsymmErrors*> gEleEffVsDb;
  std::map<int, TGraphAsymmErrors*> gEleEffVsBsChargeCorr;
  std::map<int, TGraphAsymmErrors*> gEleEffVsGenChargeCorr;

  TGraphAsymmErrors* gEleEffVsPtTest;
  
  std::vector<double> vPtBinEdges;
  std::vector<double> vMVACatBinEdges;
  
  unsigned int HZZV1IDEleBit;
  unsigned int HZZV2IDEleBit;
  unsigned int MVANILIDEleBit;
  unsigned int MVAILIDEleBit;

  // Used for histograms
  std::map<int, std::string> genEleHistMap;
  std::map<int, std::string> genEleAllMap;
  std::map<int, std::string> genEleFromBMap;
  std::map<int, std::string> genEleMatchedMap;
  std::map<int, std::string> genEleBsChargeCorrMap;

  // Used for efficiencies (TGraphAsymmErrors)
  std::map<int, std::string> genEleGraphMap;

  
  std::map<int, std::string> eleHistMap;
  std::map<int, std::string> eleAllMap;
  std::map<int, std::string> eleGenMap;
  std::map<int, std::string> elePVMap;
  std::map<int, std::string> eleSelMap;
  std::map<int, std::string> eleBsChargeCorrMap;
  unsigned int eleSelShift;

  // Used for efficiencies (TGraphAsymmErrors)
  std::map<int, std::string> eleGraphMap;
  
  // dummy copy constructor and assignment
  ElectronEfficiencyAnalyzer           ( const ElectronEfficiencyAnalyzer& );
  ElectronEfficiencyAnalyzer& operator=( const ElectronEfficiencyAnalyzer& );

};


#endif

