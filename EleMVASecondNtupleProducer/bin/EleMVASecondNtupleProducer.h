#ifndef EleMVASecondNtupleProducer_H
#define EleMVASecondNtupleProducer_H
 
#include "TH1.h"
#include "TH2.h"
#include "PDAnalysis/Ntu/interface/MGBaseAnalyzer.h"
#include "PDAnalysis/Ntu/interface/MGSelector.h"
#include "PDAnalysis/Ntu/interface/MGGenTools.h"
#include "MGTools/PlottingTools/interface/PlottingTools.h"

class EleMVASecondNtupleWriter;

class EleMVASecondNtupleProducer: 
    public MGBaseAnalyzer,
    public virtual MGSelector,
    public virtual MGGenTools, 
    public virtual PlottingTools 
{
 public:

  EleMVASecondNtupleProducer() : MGBaseAnalyzer("EleMVASecondNtupleProducer") {};
  virtual ~EleMVASecondNtupleProducer();

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

  EleMVASecondNtupleWriter* tWriter;
  
  std::string secondNtupleFileName;
  std::string tightSelection;
  
  float trueBMass;
  
  int nConeIterations;

 private:
  
  TH2D* hEleBTrkDistance;
  TCanvas* cEleBTrkDistance;
  
  TH1D* hWeightedBMass;
  TCanvas* cWeightedBMass;  
   
  TH1D* hEleConeDistance;
  TCanvas* cEleConeDistance;
  
  // dummy copy constructor and assignment
  EleMVASecondNtupleProducer           ( const EleMVASecondNtupleProducer& );
  EleMVASecondNtupleProducer& operator=( const EleMVASecondNtupleProducer& );

};

#endif
