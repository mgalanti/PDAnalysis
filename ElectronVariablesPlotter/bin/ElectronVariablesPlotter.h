#ifndef ElectronVariablesPlotter_H
#define ElectronVariablesPlotter_H
 
#include "TH1.h"
#include "TH2.h"
#include "PDAnalysis/Ntu/interface/MGBaseAnalyzer.h"
#include "PDAnalysis/Ntu/interface/MGSelector.h"
#include "PDAnalysis/Ntu/interface/MGGenTools.h"
#include "MGTools/PlottingTools/interface/PlottingTools.h"

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
   
   std::vector<TH1D*> vhEleIDTypeDecision;
   std::vector<TCanvas*> vcEleIDTypeDecision;
   
   std::vector<TH1D*> vhEleMVAType;
   std::vector<TCanvas*> vcEleMVAType;

   std::vector<TH1D*> vhEleUserFloat;
   std::vector<TCanvas*> vcEleUserFloat;
   
   std::vector<TH1D*> vhEleUserInt;
   std::vector<TCanvas*> vcEleUserInt;
   
   std::map<int,std::string> eleMVATypeMap;
   std::map<int,std::string> eleUserFloatMap;
   std::map<int,std::string> eleUserIntMap;
   
  
   
   
  // dummy copy constructor and assignment
  ElectronVariablesPlotter           ( const ElectronVariablesPlotter& );
  ElectronVariablesPlotter& operator=( const ElectronVariablesPlotter& );

};

#endif
