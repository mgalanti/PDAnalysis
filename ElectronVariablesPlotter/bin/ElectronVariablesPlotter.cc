#include <iostream>
#include <sstream>
#include <string>
#include <bitset>
#include <math.h>

#include "ElectronVariablesPlotter.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCollection.h"
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
  
  processName = selectionSubStrings[0];
  
  chargeCorrCutName = "";
  chargeCorrCut = 0;
  if(selectionSubStrings.size() > 2)
  {
    chargeCorrCutName = selectionSubStrings[2];
    if(chargeCorrCutName.compare("pos2") == 0)
      chargeCorrCut = 2;
    if(chargeCorrCutName.compare("pos") == 0)
      chargeCorrCut = 1;
    if(chargeCorrCutName.compare("neg") == 0)
      chargeCorrCut = -1;
    if(chargeCorrCutName.compare("neg2") == 0)
      chargeCorrCut = -2;
  }
  
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
  
  autoSavedObject = hEleBChargeCorr = Create1DHistogram<TH1D>("hEleBChargeCorr", "Charge correlation between e and B", 5, -2.5, 2.5, "Correlation", "N. electrons");
  
  int counter = 0;
  for(auto itEleMVAType : eleMVATypeMap)
  {
    std::string name = itEleMVAType.second;
    std::string histoName = "h" + name;
    int nBins = 100;
    int xMin = -1.;
    int xMax = 1.;
    if(histoName.find("RawValues") != std::string::npos)
    {
      nBins = 100;
      xMin = -5.;
      xMax = 5.;
    }
    TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), name.c_str(), nBins, xMin, xMax, histoName.c_str(), "N. electrons");
    vhEleMVAType.push_back(histo);
    autoSavedObject = vhEleMVAType[counter++];
  }
  
  counter = 0;
  for(auto itEleUserFloat : eleUserFloatMap)
  {
    std::string name = itEleUserFloat.second;
    std::string histoName = "h" + name;
    TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), name.c_str(), 100, 0., 100., histoName.c_str(), "N. electrons");
    vhEleUserFloat.push_back(histo);
    autoSavedObject = vhEleUserFloat[counter++];
  }

  counter = 0;
  for(auto itEleUserInt : eleUserIntMap)
  {
    std::string name = itEleUserInt.second;
    std::string histoName = "h" + name;
    int nBins = 8;
    double xMin = -0.5;
    double xMax = 7.5;
    if(histoName.find("cutBased") != std::string::npos)
    {
      nBins = 4097;
      xMin = -0.5;
      xMax = 4096.5;
    }
    TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), name.c_str(), nBins, xMin, xMax, histoName.c_str(), "N. electrons");
    vhEleUserInt.push_back(histo);
    autoSavedObject = vhEleUserInt[counter++];
  }
  
  for(auto eleFloat : mAllEleFloats)
  {
    std::string name = eleFloat.first;
    std::string histoName  = "h" + name;
    int nBins = 100;
    float xMin = 0.;
    float xMax = 100.;
    if(name.find("Eta") != std::string::npos)
    {
      xMax = 3.;
      if(name.find("Abs") != std::string::npos)
      {
        xMin = 0.;
      }
      else
      {
        xMin = -3.;
      }
    }
    else if(name.find("Phi") != std::string::npos)
    {
      xMin = -constants::PI;
      xMax = constants::PI;
    }
    else if(name.find("Dxy") != std::string::npos)
    {
      xMin = 0;
      xMax = 2.;
    }
    else if(name.find("Dz") != std::string::npos)
    {
      xMin = 0;
      xMax = 20.;
    }
    else if(name.find("Db") != std::string::npos)
    {
      xMin = 0;
      xMax = 10.;
    }
    else if(name.find("Exy") != std::string::npos)
    {
      xMin = 0.;
      xMax = 1.;
    }
    else if(name.find("Ez") != std::string::npos)
    {
      xMin = 0.;
      xMax = 1.;
    }
    else if(name.find("Conv") != std::string::npos)
    {
      xMin = 0;
      xMax = 20.;
    }
    else if(name.find("CtfGsfOv") != std::string::npos)
    {
      xMin = 0.;
      xMax = 1.;
    }
    else if(name.find("FBrem") != std::string::npos)
    {
      xMin = 0.;
      xMax = 1.;      
    }

    if(name.find("F5x5S") != std::string::npos)
    {
      xMin = 0.;
      xMax = 0.1;
    }
    if(name.find("Width") != std::string::npos)
    {
      xMin = 0.;
      xMax = 0.5;
    }
    if(name.find("Prob") != std::string::npos)
    {
      xMin = 0.;
      xMax = 1.;
    }

    
    TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), name.c_str(), nBins, xMin, xMax, name.c_str(), "N. electrons");
    mhEleVariables.insert(std::make_pair(name,histo));
    autoSavedObject = mhEleVariables[name];
  }
  
  for(auto eleInt : mAllEleInts)
  {
    std::string name = eleInt.first;
    std::string histoName  = "h" + name;
    int nBins = 100;
    float xMin = -0.5;
    float xMax = 99.5;
    if(name.find("Charge") != std::string::npos)
    {
      nBins = 3;
      xMin = -1.5;
      xMax = 1.5;
    }
    else if(name.find("Pattern") != std::string::npos)
    {
      nBins = 200;
      xMin = 0.;
      xMax = 65535;
    }
    else if(name.find("ID") != std::string::npos)
    {
      nBins = 200;
      xMin = 0;
      xMax = 1e9;
    }
    else if(name.find("Props") != std::string::npos)
    {
      nBins = 200;
      xMin = 0.;
      xMax = 1048576.;
    }
    
    TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), name.c_str(), nBins, xMin, xMax, name.c_str(), "N. electrons");
    mhEleVariables.insert(std::make_pair(name,histo));
    autoSavedObject = mhEleVariables[name];
  }
  
  for(auto eleBool : mAllEleBools)
  {
    std::string name = eleBool.first;
    std::string histoName  = "h" + name;
    int nBins = 2;
    float xMin = -0.5;
    float xMax = 1.5;
    TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), name.c_str(), nBins, xMin, xMax, name.c_str(), "N. electrons");
    mhEleVariables.insert(std::make_pair(name,histo));
    autoSavedObject = mhEleVariables[name];
  }
  
  for(auto nObject : mAllNObjects)
  {
    std::string name = nObject.first;
    if(name.compare("nElectrons") == 0)
    {
      std::string histoName  = "h" + name;
      int nBins = 21;
      float xMin = -0.5;
      float xMax = 20.5;
      TH1D* histo = Create1DHistogram<TH1D>(histoName.c_str(), name.c_str(), nBins, xMin, xMax, name.c_str(), "N. electrons");
      mhEleVariables.insert(std::make_pair(name,histo));
      autoSavedObject = mhEleVariables[name];
      
      break;
    }
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
  
  if(entry == 0)
  {
    // For BsToJPsiPhi we require generator information to work
    if(processName.compare("BsToJPsiPhi") == 0 && (!use_gen || !has_gen))
    {
      std::cout << "ElectronVariablesPlotter::analyze(): E R R O R ! processName is " << processName << " but generator information";
      if(!has_gen)
        std::cout << " is not present in the ntuple";
      if(!use_gen && !has_gen)
        std::cout << " and";
      if(!use_gen)
        std::cout << " is not used";
      std::cout << "!\n";
      std::cout << "                                     Please fix the configuration parameters.\n";
      std::cout << "                                     Exiting...\n";
      exit(1);
    }
  }
  
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
  
  // Select the event
  bool evtSelected = SelectEvent();
  
  if(!evtSelected)
    return false;
 
  // Do something with the event and the selected objects...
  
  // Retrieve all the selected objects and check that everything is ok
  int iSelObject = 0;
  int iBestPV = -1;
  int iBestB = -1;
  int iBestEle = -1;
  for (auto itSelObjects = selectedObjects.begin(); itSelObjects != selectedObjects.end(); itSelObjects++)
  {
    //std::cout << "Event selection: selected object #" << iSelObject << ": type = " << itSelObjects->first << ", index = " << itSelObjects->second << std::endl;
    if(itSelObjects->first == PDEnumString::recPV)
      iBestPV = itSelObjects->second;
    if(itSelObjects->first == PDEnumString::recSvt)
      iBestB = itSelObjects->second;
    if(itSelObjects->first == PDEnumString::recElectron)
      iBestEle = itSelObjects->second;
    iSelObject++;
  }
    
  // Depending on the selection string provided in the configuration, it is not granted 
  // that we have all the needed objects at this point, even if the event passed the selection.
  // Thus, let's check explicitly.  
  if(iBestPV == -1 || iBestB == -1 || iBestEle == -1)
  {
    std::cout << "E R R O R ! Event passed the selection \"" << evtSelection << "\" but not all the needed objects are available!\n";
    std::cout << "            iBestPV = " << iBestPV << ", iBestB = " << iBestB << ", iBestEle = " << iBestEle << std::endl;
    std::cout << "            Please fix the configuration file or the MGSelector::SelectEvent() code to avoid this error.";
    std::cout << "            Exiting...\n";
    exit(1);
  }
  
  // Check charge correlation (this depends on the processName and on whether the gen information is used or not)
  std::vector<int> tracksFromB = tracksFromSV(iBestB);
  int chargeB = 0;
  for(auto iTrack : tracksFromB)
  {
    chargeB += trkCharge->at(iTrack);
  }
  int chargeEle = eleCharge->at(iBestEle);
  
  // If we don't use generator information, then that's it
  int chargeCorr = chargeEle * chargeB;
  
  // Otherwise let's retrieve the charge correlation from gen
  if(use_gen)
  {
    int iMatchedGenEle = GetClosestGen(elePt->at(iBestEle), eleEta->at(iBestEle), elePhi->at(iBestEle));
//     int idMatchedGenEle = 0;
//     if (iMatchedGenEle >= 0)
//       idMatchedGenEle = genId->at(iMatchedGenEle);
    int iGenB = GetClosestGenNoLL(iBestB);
//     int idGenB = 0;
//     if(iGenB > 0)
//     {
//       idGenB = genId->at(iGenB);
//     }
    if(processName.compare("BsToJPsiPhi") == 0)
    {
      if(iMatchedGenEle)
        chargeCorr = GetGenLepBsChargeCorrelation(iMatchedGenEle, iGenB);
      if(!chargeCorr)
      {
        chargeCorr = GetBsChargeCorrelation(chargeEle, iGenB);
        // In this case, correlation value goes to the +/-2 bins
        chargeCorr*=2;
      }
    }
    else if(processName.compare("BuToJPsiK") == 0)
    {
      int recoChargeCorr = chargeCorr;
      chargeCorr = GetGenLepBuChargeCorrelation(iMatchedGenEle, iGenB);
      if(!chargeCorr)
      {
        chargeCorr = GetBuChargeCorrelation(chargeEle, iGenB);
        // In this case, correlation value goes to the +/-2 bins
        chargeCorr*=2;
      }
      if(recoChargeCorr * chargeCorr < 0)
      {
        std::cout << "W A R N I N G ! Bu-ele charge correlation gives different sign if retrieved from gen or from reco!\n";
      }
    }
  }
  
  // Now select based on the charge correlation
  if(chargeCorrCut && chargeCorr != chargeCorrCut)
    return false;

  hEleBChargeCorr->Fill(chargeCorr);  
    
  mhEleVariables["nElectrons"]->Fill(*(mAllNObjects["nElectrons"]));
  
  for (int iUserInfo = 0; iUserInfo < nUserInfo; iUserInfo++)
  {
    if(useObjType->at(iUserInfo) == PDEnumString::recElectron && useObjIndex->at(iUserInfo) == iBestEle)
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
          if(nameEleMVA.compare(histo->GetTitle()) == 0)
            histo->Fill(infoValue);
        }
      }
      auto eleUserFloatMatch = eleUserFloatMap.find(infoType);
      if(eleUserFloatMatch != eleUserFloatMap.end())
      {
        std::string nameEleUserFloat = eleUserFloatMatch->second;
        for(auto histo : vhEleUserFloat)
        {
          if(nameEleUserFloat.compare(histo->GetTitle()) == 0)
            histo->Fill(infoValue);
        }
      }
      auto eleUserIntMatch = eleUserIntMap.find(infoType);
      if(eleUserIntMatch != eleUserIntMap.end())
      {
        std::string nameEleUserInt = eleUserIntMatch->second;
        for(auto histo : vhEleUserInt)
        {
          if(nameEleUserInt.compare(histo->GetTitle()) == 0)
            histo->Fill(infoValue);
        }
      }
    }
  }

  // Now fill the electron variables
  for(auto eleFloat : mAllEleFloats)
  {
    std::string name = eleFloat.first;
    std::vector<float>** value = eleFloat.second;
    //       std::cout << "Electron # " << iBestEle << ": variable with name " << name << " has value " << (*value)->at(iElectron) << std::endl;
    mhEleVariables[name]->Fill((*value)->at(iBestEle));
  }
  for(auto eleInt : mAllEleInts)
  {
    std::string name = eleInt.first;
    std::vector<int>** value = eleInt.second;
    //       std::cout << "Electron # " << iBestEle << ": variable with name " << name << " has value " << (*value)->at(iElectron) << std::endl;
    mhEleVariables[name]->Fill((*value)->at(iBestEle));
  }
  for(auto eleBool : mAllEleBools)
  {
    std::string name = eleBool.first;
    std::vector<bool>** value = eleBool.second;
    //       std::cout << "Electron # " << iBestEle << ": variable with name " << name << " has value " << (*value)->at(iElectron) << std::endl;
    mhEleVariables[name]->Fill((*value)->at(iBestEle));
  }
  
  return true;
}


void ElectronVariablesPlotter::endJob() 
{
  // This runs after the event loop
  
  autoSavedObject = cEleBChargeCorr = CreateCanvas("cEleBChargeCorr", 0, 21, 1, false, false, hEleBChargeCorr);
  
  TCanvas* cTemp;
  
  for(auto histo : JoinMany<TH1D*>({vhEleMVAType, vhEleUserFloat, vhEleUserInt}))
  {
    std::string name = histo->GetName();
    name = name.substr(1);
    std::string canvasName = "c" + name;
    std::cout << "ElectronVariablesPlotter::endJob(): I N F O. Creating TCanvas " << canvasName << std::endl;
    autoSavedObject = cTemp = CreateCanvas(canvasName.c_str(), 0, 21, 1, false, false, histo);
  }

  for(auto hEleVariable : mhEleVariables)
  {
    std::string name = hEleVariable.first;
    std::string canvasName = "c" + name;
    std::cout << "ElectronVariablesPlotter::endJob(): I N F O. Creating TCanvas " << canvasName << std::endl;
    autoSavedObject = cTemp = CreateCanvas(canvasName.c_str(), 0, 21, 1, false, false, hEleVariable.second);
  }
    
  return;
}
