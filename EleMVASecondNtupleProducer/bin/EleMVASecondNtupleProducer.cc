#include <iostream>
#include <sstream>
#include <string>
#include <bitset>
#include <math.h>

#include "EleMVASecondNtupleProducer.h"

#include "PDAnalysis/Ntu/interface/PDSecondNtupleWriter.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"

using namespace std;



EleMVASecondNtupleProducer::~EleMVASecondNtupleProducer() {
}



void EleMVASecondNtupleProducer::beginJob() {

  MGBaseAnalyzer::beginJob();


  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  getUserParameter( "verbose", verbose );

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor
  
  return;
}


void EleMVASecondNtupleProducer::book() 
{

  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility

  return;
}


void EleMVASecondNtupleProducer::reset() {
// automatic reset
  autoReset();
  return;
}


bool EleMVASecondNtupleProducer::analyze( int entry, int event_file, int event_tot ) {

  if ( verbose ) {
    cout << " +++++++++++++++++++++++++++ " << endl;
    cout << "entry: "
         << entry << " " << event_file << " " << event_tot << endl;
    cout << "run: " <<   runNumber << " , "
         << "evt: " << eventNumber << endl;
  }
  else {
//    if ( !( event_file % 10000 ) || !( event_tot % 10000 ) )
    if ( !( event_tot % 10000 ) && event_tot )
      cout << event_file << " " << event_tot << endl;
  }

  // Main per-event analysis code goes here
  
  // Event and object selection code
  // Method SelectEvent(...) is defined in class MGSelector
  //    evtSelection string is set in the configuration file
  //    selectedObjects is filled by SelectEvent to store pointers to all selected objects in the event
  //       The format is a pair<int, int>, 
  //       where the first integer is the object type 
  //       (as defined in the PDEnumString::recoObject enum)
  //       and the second one is the object index
  std::vector<std::pair<int,int> > selectedObjects;
  
  bool evtSelected = SelectEvent(evtSelection.c_str(), selectedObjects);
  
  if(!evtSelected)
    return false;
 
  // Do something with the event and the selected objects...
    
  return true;
}


void EleMVASecondNtupleProducer::endJob() 
{
  // This runs after the event loop
  
  return;
}
