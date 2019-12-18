#ifndef EleMVACalibrationProducer_h
#define EleMVACalibrationProducer_h

#include <string>
#include <iostream>

#include "MGTools/PlottingTools/interface/PlottingTools.h"
#include "PDAnalysis/EleMVASecondNtupleProducer/interface/EleMVASecondNtupleData.h"
#include "NtuTool/Read/interface/TreeReader.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TMVA/Reader.h"

class EleMVACalibrationProducer: 
  public EleMVASecondNtupleData, 
  public TreeReader,
  public PlottingTools
{

 public:

  EleMVACalibrationProducer() {}
  virtual ~EleMVACalibrationProducer() {}

  virtual void beginJob();
  
//   virtual void setupNtuple();
  
//   virtual bool getEntry(int ientry);

  virtual void book();
  
  virtual bool analyze( int entry, int event_file, int event_tot );
  
  virtual void endJob();
  
  std::pair<double, double> CountEventsWithFit(TH1 *hist);

 private:
   
  bool verbose;

  std::string sampleName;
  
  TH1D* hNGenB;
  
  TH2D* hIsTagRightVsTagTruth;
  
  std::vector<TH1D*>* vhMassCalRightTag;
  std::vector<TH1D*>* vhMassCalWrongTag;
  
  TH1D* hMva;
  TH1D* hMvaRightTag;
  TH1D* hMvaWrongTag;
  TH1D* hMassTot;
  TH1D* hMassAllTag;
  TH1D* hMassRightTag;
  TH1D* hMassWrongTag;
  TH1D* hMassNoTag;
  
  TGraph* gNRightTagVsEleIdCut;
  TGraph* gNWrongTagVsEleIdCut;
  TGraph* gNNoTagVsEleIdCut;
  TGraph* gBaseEffVsEleIdCut;
  TGraph* gBaseRightTagEffVsEleIdCut;
  TGraph* gBaseWrongTagEffVsEleIdCut;
  TGraph* gBaseMistagVsEleIdCut;
  TGraph* gBaseDilutionVsEleIdCut;
  TGraph* gBasePowerVsEleIdCut;
  
  TGraph* gNRightTagVsBaseEff;
  TGraph* gNWrongTagVsBaseEff;
  TGraph* gNNoTagVsBaseEff;
  TGraph* gBaseEffVsBaseEff;
  TGraph* gBaseRightTagEffVsBaseEff;
  TGraph* gBaseWrongTagEffVsBaseEff;
  TGraph* gBaseMistagVsBaseEff;
  TGraph* gBaseDilutionVsBaseEff;
  TGraph* gBasePowerVsBaseEff;

  TGraph* gBaseWrongTagEffVsBaseRightTagEff;

  std::vector<double> vEleIdcuts;
  std::vector<double> vEleIdcutsLog;
  
  std::vector<TH1D*> vhMassRightTag;
  std::vector<TH1D*> vhMassWrongTag;

  // These are for the baseline cuts used in the calibration job
  TH1D* hElePt;
  TH1D* hEleEta;
  TH1D* hEleDxy;
  TH1D* hEleExy;
  TH1D* hEleDz;
  TH1D* hEleEz;
  TH1D* hEleIDNIV2Val;
  TH1D* hEleDRB;
  TH1D* hElePFIsoScaled;
  TH1D* hEleConeCleanPt;
  TH1D* hEleConeCleanPtRel;
  TH1D* hEleConeCleanDR;
  TH1D* hEleConeCleanEnergyRatio;
  TH1D* hEleConeCleanQ;
// FIXME: These are not implemented yet
//   // These are for all the EleId cuts defined in vEleIdcuts;
//   std::vector<TH1D*> vhElePtVsEleIdCut;
//   std::vector<TH1D*> vhEleEtaVsEleIdCut;
//   std::vector<TH1D*> vhEleDxyVsEleIdCut;
//   std::vector<TH1D*> vhEleExyVsEleIdCut;
//   std::vector<TH1D*> vhEleDzVsEleIdCut;
//   std::vector<TH1D*> vhEleEzVsEleIdCut;
//   std::vector<TH1D*> vhEleIDNIV2ValVsEleIdCut;
//   std::vector<TH1D*> vhEleDRBVsEleIdCut;
//   std::vector<TH1D*> vhElePFIsoScaledVsEleIdCut;
//   std::vector<TH1D*> vhEleConeCleanPtVsEleIdCut;
//   std::vector<TH1D*> vhEleConeCleanPtRelVsEleIdCut;
//   std::vector<TH1D*> vhEleConeCleanDRVsEleIdCut;
//   std::vector<TH1D*> vhEleConeCleanEnergyRatioVsEleIdCut;
//   std::vector<TH1D*> vhEleConeCleanQVsEleIdCut;
  
  std::vector<TH1D*> vhCalParsRightTagVsMvaScore;
  std::vector<TH1D*> vhCalParsWrongTagVsMvaScore;  
  TH1D* hCalFitStatusRightTagVsMvaScore;
  TH1D* hCalFitStatusWrongTagVsMvaScore;
  TH1D* hCalFitChi2NDOFRightTagVsMvaScore;
  TH1D* hCalFitChi2NDOFWrongTagVsMvaScore;
  
  std::string treeListName;
  std::string inputTreeName;
  std::string mvaInputPath;
  std::string mvaMethod;
  std::string process;
  std::string dirPath;
  bool useTightSelection;
  bool weightControlPlots;
  bool isData = false;
  bool useSyst;
  bool fixToTot;
  bool fixToAllTag;
  bool computeLinearitySyst;
  bool linearitySystSecOrd;
  bool linearitySystThirdOrd;
  int nBinsCal;
  int nBinsMva;
  double minMass, maxMass, x1MassSB, x2MassSB;
  double meanTotMass, sigma1TotMass, sigma2TotMass, fracTotMass;
  double erfcWidthTotMass, erfcShiftTotMass;
  double meanAllTagMass, sigma1AllTagMass, sigma2AllTagMass, fracAllTagMass;
  double erfcWidthAllTagMass, erfcShiftAllTagMass;
  int nMassBins = 50;
  bool writeOutput = true;
  double systematics2017 = 0.;
  double systematics2018 = 0.;
  double systematics = 0.;
  double eleIdWP;
  double elePtWP;
  double eleDzWP;
  double eleDxyWP;
  double eleDRBWP;
  double mvaValueMinCutOff;
  double mvaValueMaxCutOff;
  
  TMVA::Reader* reader;
  
  std::vector<double> evtW; // per-event mistag rate
  double totalP;  // total tagging power
  double totalPBinned;
  double totalPBinnedErr;
  double binSizeCal;
  
  std::vector<double> vWCalc;
  std::vector<double> vWCalcLowEdge;
  std::vector<double> vWCalcHighEdge;
  
  TFitResultPtr fitResult; // Stores the result of the last fit performed
  
  
  EleMVACalibrationProducer           ( const EleMVACalibrationProducer& a );
  EleMVACalibrationProducer& operator=( const EleMVACalibrationProducer& a );
};

#endif
