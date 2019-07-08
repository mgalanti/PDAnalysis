#include <iostream>
#include <sstream>
#include <string>
#include <bitset>
#include <math.h>

#include "OSElectronFlavorTagger.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"

using namespace std;



OSElectronFlavorTagger::~OSElectronFlavorTagger() {
}



void OSElectronFlavorTagger::beginJob() {

  MGBaseAnalyzer::beginJob();
  MGSelector::SetSelectionString(evtSelection);
  
  initializeOsElectronMvaReader();


  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor
  
  return;
}


void OSElectronFlavorTagger::book() 
{

  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility

  return;
}


void OSElectronFlavorTagger::reset() {
// automatic reset
  autoReset();
  return;
}


bool OSElectronFlavorTagger::analyze( int entry, int event_file, int event_tot ) {

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
  if(entry == 0)
    checkBranches();
  
  // Main per-event analysis code goes here
  
  // Event and object selection code
  // Method SelectEvent() is defined in class MGSelector
  // It depends on two variables:
  //    evtSelection is a string set in the configuration file, which determines the selection algorithm and parameters
  //    selectedObjects is declared in MGSelector.h. It is filled by SelectEvent to store pointers to all selected objects in the event
  //       The format of selectedObjects is a vector<pair<int, int> >, 
  //       where the first integer is the object type (as defined in the PDEnumString::recoObject enum)
  //       and the second one is the object index
  bool evtSelected = SelectEvent();
  // In alternative, one can use custom selection and vector of selected objects:
  // std::string mySelection = "mySelectionString";  
  // std::vector<std::pair<int,int> > mySelectedObjects;
  // bool myEvtSelected = SelectEvent(mySelection.c_str(), mySelectedObjects);
  
  if(!evtSelected)
    return false;
  
  int iSelObject = 0;
  int iBestPV = -1;
  int iBestB = -1;
  int iBestEle = -1;
  // Do something with the event and the selected objects...
  for (auto itSelObjects = selectedObjects.begin(); itSelObjects != selectedObjects.end(); itSelObjects++)
  {
    //     std::cout << "Event selection: selected object #" << iSelObject << ": type = " << itSelObjects->first << ", index = " << itSelObjects->second << std::endl;
    if(itSelObjects->first == PDEnumString::recPV)
      iBestPV = itSelObjects->second;
    if(itSelObjects->first == PDEnumString::recSvt)
      iBestB = itSelObjects->second;
    if(itSelObjects->first == PDEnumString::recElectron)
      iBestEle = itSelObjects->second;
    iSelObject++;
  }
  
  setObjectIndexes(iBestEle, iBestB, iBestPV);
  bool tagOk = makeOsElectronTagging();
  std::cout << "tagOk = " << tagOk << std::endl;
  
  std::cout << "osElectronTagDecision = " << getOsElectronTag() << std::endl;
  std::cout << "osElectronTagMvaValue = " << getOsElectronTagMvaValue() << std::endl;
  
  float fixedGridRhoFastJetAll = GetFixedGridRhoFastJetAll();
  float fixedGridRhoFastJetAllCalo = GetFixedGridRhoFastJetAllCalo();
  
  std::cout << "fixedGridRhoFastJetAll = " << fixedGridRhoFastJetAll << std::endl;
  std::cout << "fixedGridRhoFastJetAllCalo = " << fixedGridRhoFastJetAllCalo << std::endl;
    
  return true;
}


void OSElectronFlavorTagger::endJob() 
{
  // This runs after the event loop
  
  return;
}
