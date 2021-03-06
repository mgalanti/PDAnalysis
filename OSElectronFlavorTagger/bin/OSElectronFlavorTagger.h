#ifndef OSElectronFlavorTagger_H
#define OSElectronFlavorTagger_H
 
#include "TH1.h"
#include "TH2.h"
#include "PDAnalysis/Ntu/interface/MGBaseAnalyzer.h"
#include "PDAnalysis/Ntu/interface/MGSelector.h"
#include "PDAnalysis/Ntu/interface/MGGenTools.h"
#include "MGTools/PlottingTools/interface/PlottingTools.h"
#include "PDAnalysis/OSElectronFlavorTagger/interface/OsElectronMvaTag.h"

class OSElectronFlavorTagger: 
    public MGBaseAnalyzer,
    public virtual MGSelector,
    public virtual MGGenTools, 
    public virtual PlottingTools,
    public virtual OsElectronMvaTag
{
 public:

  OSElectronFlavorTagger() : MGBaseAnalyzer("OSElectronFlavorTagger") {};
  virtual ~OSElectronFlavorTagger();

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
  
  int nMisTags;
  int nAllTags;
  int nAllBEvts;
  float sumMisTagProb;

 private:
  
  // dummy copy constructor and assignment
  OSElectronFlavorTagger           ( const OSElectronFlavorTagger& );
  OSElectronFlavorTagger& operator=( const OSElectronFlavorTagger& );

};

#endif
