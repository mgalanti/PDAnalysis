#include <iostream>
#include <sstream>
#include <string>
#include <bitset>
#include <math.h>

#include "TSystem.h"

#include "EleMVASecondNtupleProducer.h"

#include "PDAnalysis/EleMVASecondNtupleProducer/interface/EleMVASecondNtupleWriter.h"

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
  MGSelector::SetSelectionString(evtSelection);


  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  getUserParameter( "verbose", verbose );

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor
  
  tWriter = new EleMVASecondNtupleWriter; // second ntuple
  
  std::string secondNtupleBaseName = getUserParameter("secondNtupleBaseName");
  
  // Set the second ntuple file name
  int tries = 0;
  secondNtupleFileName = secondNtupleBaseName + "__" + sampleName + "__" + evtSelection + "__" + std::to_string(tries++) + ".root";
  std::cout << "EleMVASecondNtupleProducer::beginJob(): setting variables:\n";
  std::cout << "                                        secondNtupleFileName = " << secondNtupleFileName << std::endl;
  bool secondNtupleOrigNameNotOk = false;
  while(!gSystem->AccessPathName(gSystem->ExpandPathName(secondNtupleFileName.c_str())))
  {
    secondNtupleOrigNameNotOk = true;
    if(tries == 1)
    {
      std::cout << "W A R N I N G! Second ntuple file \"" << secondNtupleFileName << "\" already exists!\n";
      std::cout << "               Trying a different name...\n";
    }
    secondNtupleFileName = secondNtupleBaseName + "__" + sampleName + "__" + evtSelection + "__" + std::to_string(tries++) + ".root";
    if(tries > 9999)
      break;
  }
  if(!gSystem->AccessPathName(gSystem->ExpandPathName(secondNtupleFileName.c_str())))
  {
    std::cout << "E R R O R! Output file \"" << secondNtupleFileName << "\" already exists!\n";
    std::cout << "           Exiting...\n";
    exit(1);
  }
  if(secondNtupleOrigNameNotOk)
    std::cout << "                                        New secondNtupleFileName = " << secondNtupleFileName << std::endl;
  
  // Actually open the second ntuple file for writing
  tWriter->open(secondNtupleFileName, "RECREATE");

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
  
  tWriter->Reset();
  
  // Event and object selection code
  // Method SelectEvent(...) is defined in class MGSelector
  //    evtSelection string is set in the configuration file
  //    selectedObjects is filled by SelectEvent to store pointers to all selected objects in the event
  //       The format is a pair<int, int>, 
  //       where the first integer is the object type 
  //       (as defined in the PDEnumString::recoObject enum)
  //       and the second one is the object index
//   std::vector<std::pair<int,int> > selectedObjects;
  
//   bool evtSelected = SelectEvent(evtSelection.c_str(), selectedObjects);
  
  bool evtSelected = SelectEvent();
  
  if(!evtSelected)
    return false;
  
 
  // Do something with the event and the selected objects...
  int iSelObject = 0;
  int iBestPV = -1;
  int iBestBs = -1;
  int iBestEle = -1;
  for (auto itSelObjects = selectedObjects.begin(); itSelObjects != selectedObjects.end(); itSelObjects++)
  {
//     std::cout << "Event selection: selected object #" << iSelObject << ": type = " << itSelObjects->first << ", index = " << itSelObjects->second << std::endl;
    if(itSelObjects->first == PDEnumString::recPV)
      iBestPV = itSelObjects->second;
    if(itSelObjects->first == PDEnumString::recSvt)
      iBestBs = itSelObjects->second;
    if(itSelObjects->first == PDEnumString::recElectron)
      iBestEle = itSelObjects->second;
    iSelObject++;
  }
  
//   int iElectron = SelectOSElectron(iBestPV, iBestBs);

  // Depending on the selection string provided in the configuration, it is not granted 
  // that we have all the needed objects at this point, even if the event passed the selection.
  // Thus, let's check explicitly.  
  if(iBestPV == -1 || iBestBs == -1 || iBestEle == -1)
  {
    std::cout << "E R R O R ! Event passed the selection \"" << evtSelection << "\" but not all the needed objects are available!\n";
    std::cout << "            iBestPV = " << iBestPV << ", iBestBs = " << iBestBs << ", iBestEle = " << iBestEle << std::endl;
    std::cout << "            Please fix the configuration file or the MGSelector::SelectEvent() code to avoid this error.";
    std::cout << "            Exiting...\n";
    exit(1);
  }
  
  (tWriter->elePt) = elePt->at(iBestEle);
  (tWriter->eleEta) = eleEta->at(iBestEle);
  (tWriter->elePhi) = elePhi->at(iBestEle);
  
  tWriter->fill();
    
  return true;
}



void EleMVASecondNtupleProducer::endJob() 
{
  // This runs after the event loop
  tWriter->close();
  
  return;
}
