#include <iostream>
#include <sstream>
#include <string>
#include <bitset>
#include <math.h>

#include "ElectronVariablesPlotter.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"

using namespace std;



ElectronVariablesPlotter::~ElectronVariablesPlotter() {
}



void ElectronVariablesPlotter::beginJob() {

  MGBaseAnalyzer::beginJob();
  MGSelector::SetSelectionString(evtSelection);
  
  eleMVATypeMap = PDEnumString::eleMVATypeMap();
  eleUserFloatMap = PDEnumString::eleUserFloatMap();
  eleUserIntMap = PDEnumString::eleUserIntMap();

  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor
  
  return;
}


void ElectronVariablesPlotter::book() 
{

  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility
  
//   for(auto itEleID : eleIDTypeMapIS)
//   {
//     std::string histoName = itEleID.second;
//     TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), histoName.c_str(), 2, -0.5, 1.5, "ID decision", "N. electrons");
//     vhEleIDTypeDecision.push_back(histo);
//   }
//   autoSavedObject = vhEleIDTypeDecision;
  
  int counter = 0;
  for(auto itEleMVAType : eleMVATypeMap)
  {
    std::string histoName = itEleMVAType.second;
    int nBins = 100;
    int xMin = -1.;
    int xMax = 1.;
    if(histoName.substr(histoName.length()-9).compare("RawValues") == 0)
    {
      nBins = 100;
      xMin = -5.;
      xMax = 5.;
    }
    TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), histoName.c_str(), nBins, xMin, xMax, histoName.c_str(), "N. electrons");
    vhEleMVAType.push_back(histo);
    autoSavedObject = vhEleMVAType[counter++];
  }
  
  counter = 0;
  for(auto itEleUserFloat : eleUserFloatMap)
  {
    std::string histoName = itEleUserFloat.second;
    TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), histoName.c_str(), 100, 0., 100., histoName.c_str(), "N. electrons");
    vhEleUserFloat.push_back(histo);
    autoSavedObject = vhEleUserFloat[counter++];
  }

  counter = 0;
  for(auto itEleUserInt : eleUserIntMap)
  {
    std::string histoName = itEleUserInt.second;
    int nBins = 8;
    double xMin = -0.5;
    double xMax = 7.5;
    if(histoName.substr(0,8).compare("cutBased") == 0)
    {
      nBins = 4097;
      xMin = -0.5;
      xMax = 4096.5;
    }
    TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), histoName.c_str(), nBins, xMin, xMax, histoName.c_str(), "N. electrons");
    vhEleUserInt.push_back(histo);
    autoSavedObject = vhEleUserInt[counter++];
  }

  return;
}


void ElectronVariablesPlotter::reset() {
// automatic reset
  autoReset();
  return;
}


bool ElectronVariablesPlotter::analyze( int entry, int event_file, int event_tot ) {

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
//   bool evtSelected = SelectEvent();
  // In alternative, one can use custom selection and vector of selected objects:
  // std::string mySelection = "mySelectionString";  
  // std::vector<std::pair<int,int> > mySelectedObjects;
  // bool myEvtSelected = SelectEvent(mySelection.c_str(), mySelectedObjects);
  
//   if(!evtSelected)
//     return false;
 
  // Do something with the event and the selected objects...
    
//   std::cout << "This event has " << useObjType->size() << " userInfo entries\n";
  for(int iElectron = 0; iElectron < nElectrons; iElectron++)
  {
//     std::cout << "Electron " << iElectron << std::endl;
    for (int iUserInfo = 0; iUserInfo < nUserInfo; iUserInfo++)
    {
      if(useObjType->at(iUserInfo) == PDEnumString::recElectron && useObjIndex->at(iUserInfo) == iElectron)
      {
//         std::cout << "   Found corresponding user info!\n";
        int infoType = useInfoType->at(iUserInfo);
        double infoValue = useInfoValue->at(iUserInfo);
//         std::cout << "      type = " << infoType << ", value = " << infoValue << std::endl;
        auto eleMVATypeMatch = eleMVATypeMap.find(infoType);
        if(eleMVATypeMatch != eleMVATypeMap.end())
        {
          std::string nameEleMVA = eleMVATypeMatch->second;
          for(auto histo : vhEleMVAType)
          {
            if(nameEleMVA.compare(histo->GetName()) == 0)
              histo->Fill(infoValue);
          }
        }
        auto eleUserFloatMatch = eleUserFloatMap.find(infoType);
        if(eleUserFloatMatch != eleUserFloatMap.end())
        {
          std::string nameEleUserFloat = eleUserFloatMatch->second;
          for(auto histo : vhEleUserFloat)
          {
            if(nameEleUserFloat.compare(histo->GetName()) == 0)
              histo->Fill(infoValue);
          }
        }
        auto eleUserIntMatch = eleUserIntMap.find(infoType);
        if(eleUserIntMatch != eleUserIntMap.end())
        {
          std::string nameEleUserInt = eleUserIntMatch->second;
          for(auto histo : vhEleUserInt)
          {
            if(nameEleUserInt.compare(histo->GetName()) == 0)
              histo->Fill(infoValue);
          }
        }
        
        
      }
    }
    
//     for(auto itEleMVAType : eleMVATypeMap)
//     {
//       int iEleMVA = itEleMVAType.first;
//       std::string nameEleMVA = itEleMVAType.second;
//       double value = userInfo(PDEnumString::recElectron, iEleMVA)[iElectron];
//       for(auto histo : vhEleMVAType)
//       {
//         if(nameEleMVA.compare(histo->GetName()) == 0)
//           histo->Fill(value);
//       }
//     }
  }
  
  return true;
}


void ElectronVariablesPlotter::endJob() 
{
  // This runs after the event loop
  
  return;
}
