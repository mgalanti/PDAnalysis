#ifndef ElectronVariablesPlotter_H
#define ElectronVariablesPlotter_H
 
#include "TH1.h"
#include "TH2.h"
#include "PDAnalysis/Ntu/interface/MGBaseAnalyzer.h"
#include "PDAnalysis/Ntu/interface/MGSelector.h"
#include "PDAnalysis/Ntu/interface/MGGenTools.h"
#include "MGTools/PlottingTools/interface/PlottingTools.h"

#include <list>

class ElectronVariablesPlotter: 
    public MGBaseAnalyzer,
    public virtual MGSelector,
    public virtual MGGenTools, 
    public virtual PlottingTools 
{
 public:

  ElectronVariablesPlotter() : MGBaseAnalyzer("ElectronVariablesPlotter") {};
  virtual ~ElectronVariablesPlotter();

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

 private:
   
   TH1D* hEleBChargeCorr;
   TCanvas* cEleBChargeCorr;
   
   std::vector<TH1D*> vhEleInfoType;
   std::vector<TCanvas*> vcEleInfoType;

   std::vector<TH1D*> vhEIdInfoType;
   std::vector<TCanvas*> vcEIdInfoType;
   
   std::map<std::string,TH1D*> mhEleVariables;
   std::map<std::string,TCanvas*> mcEleVariables;

   std::map<int,std::string> eleInfoTypeMap;
   std::map<int,std::string> eIdInfoTypeMap;
      
   std::string processName;   
   std::string chargeCorrCutName;
   int chargeCorrCut;
   
   void bookEleVariableHisto(const std::string name);
   
  // dummy copy constructor and assignment
  ElectronVariablesPlotter           ( const ElectronVariablesPlotter& );
  ElectronVariablesPlotter& operator=( const ElectronVariablesPlotter& );

};

#endif
