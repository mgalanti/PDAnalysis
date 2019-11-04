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
  std::vector<TH1D*>* vhMassCalRightTag;
  std::vector<TH1D*>* vhMassCalWrongTag;
  
  TH1D* hMva;
  TH1D* hMvaRightTag;
  TH1D* hMvaWrongTag;
  TH1D* hMassTot;
  TH1D* hMassRightTag;
  TH1D* hMassWrongTag;
  TH1D* hMassNoTag;

  std::vector<double> vEleIdcuts;
  
  std::vector<TH1D*> vhMassRightTag;
  std::vector<TH1D*> vhMassWrongTag;
  
  std::string treeListName;
  std::string inputTreeName;
  std::string mvaInputPath;
  std::string mvaMethod;
  std::string process;
  std::string dirPath;
  bool useTightSelection;
  bool isData = false;
  bool useSyst;
  int nBinsCal;
  int nBinsMva;
  double minMass, maxMass, x1MassSB, x2MassSB;
  double meanTotMass, sigma1TotMass, sigma2TotMass;
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
  
  
  EleMVACalibrationProducer           ( const EleMVACalibrationProducer& a );
  EleMVACalibrationProducer& operator=( const EleMVACalibrationProducer& a );
};

#endif
