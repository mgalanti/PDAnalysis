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
  
  nMisTags = 0;
  nAllTags = 0;
  nAllBEvts = 0;
  sumMisTagProb = 0;
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

  if(iBestPV < 0 || iBestB < 0)
    return false;
  
  nAllBEvts++;
  
  if(!evtSelected)
    return false;
  
  float eleIDNIV2Val = -1;
  for (int iUserInfo = 0; iUserInfo < nUserInfo; iUserInfo++)
  {
    if(useObjType->at(iUserInfo) == PDEnumString::recElectron && useObjIndex->at(iUserInfo) == iBestEle)
    {
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2Values)
      {
        eleIDNIV2Val = useInfoValue->at(iUserInfo);
      }
    }
  }
  
  if(eleIDNIV2Val > -1.)
    return false;
  
  setObjectIndexes(iBestEle, iBestB, iBestPV);
  bool tagOk = makeOsElectronTagging();
  std::cout << "tagOk = " << tagOk << std::endl;
  
  float tagDecision = getOsElectronTag();
  float tagMvaValue = getOsElectronTagMvaValue();
  float misTagProb = 1 - tagMvaValue;
  
  std::cout << "tagDecision = " << tagDecision << std::endl;
  std::cout << "tagMvaValue = " << tagMvaValue << std::endl;
  std::cout << "misTagProb = " << misTagProb << std::endl;
    
//   float fixedGridRhoFastJetAll = GetFixedGridRhoFastJetAll();
//   float fixedGridRhoFastJetAllCalo = GetFixedGridRhoFastJetAllCalo();
  
//   std::cout << "fixedGridRhoFastJetAll = " << fixedGridRhoFastJetAll << std::endl;
//   std::cout << "fixedGridRhoFastJetAllCalo = " << fixedGridRhoFastJetAllCalo << std::endl;
  
  std::vector<int> tracksFromB = tracksFromSV(iBestB);
  TLorentzVector pB = GetTLorentzVectorFromJPsiX(iBestB);
  int iBestJPsi = (subVtxFromSV(iBestB)).at(0);
  std::vector<int> tracksFromJPsi = tracksFromSV(iBestJPsi);
 
  // Truth information for B
  int iGenB = -1;
  int nGenB = -1;
  int idGenB = 0;
  int genBMixStatus = -1;
  int evtWeight = 1;

  // Other variables needed for the ntuple production
  std::vector<int> allBHadrons;
  std::vector<int> longLivedBHadrons;

  // Use generator information if available
  if(has_gen)
  {
    // Truth information for B
    allBHadrons = GetAllBHadrons(true);  // true removes particles with mix status = 2
    longLivedBHadrons = GetAllLongLivedBHadrons(true); // true removes particles with mix status = 2
    iGenB = GetClosestGenInList(pB.Pt(), pB.Eta(), pB.Phi(), longLivedBHadrons, 0.4, 0.4);
    nGenB = longLivedBHadrons.size();
    // Do not consider events where the signal-side is not matched to a gen B
    if(iGenB < 0)
      return false;

    idGenB = genId->at(iGenB);

    // Do not consider events where the signal B is matched to a hadron different from the one looked for
    if(selectionSubStrings[0].compare("BsToJPsiPhi") == 0 && abs(idGenB) != 531)
    {
      std::cout << "I N F O: Event rejected - looking for BsToJPsiPhi but signal B is matched to a hadron with id = " << idGenB << ".\n";
      return false;
    }
    if(selectionSubStrings[0].compare("BuToJPsiK") == 0 && abs(idGenB) != 521)
    {
      std::cout << "I N F O: Event rejected - looking for BuToJPsiK but signal B is matched to a hadron with id = " << idGenB << ".\n";
      return false;
    }
    // Check whether the matched B has mixed or will mix
    genBMixStatus = GetMixStatus(iGenB);
    
    if(genBMixStatus == 2)
    {
      std::cout << "I N F O: Event rejected - mix status of generated B is 2.\n";
      return false;
    }
    
    // If the matched B hadron comes from mixing, then go back to the pre-mixing one
    if(genBMixStatus == 1)
    {
      iGenB = RecursiveLookForMotherIds(iGenB, {-idGenB});
      idGenB = genId->at(iGenB);
    }

    // Change event weight according to how many "interesting" B hadrons are there
    for(auto iGen: longLivedBHadrons)
    {
      if(iGen == iGenB)
        continue;
      if(abs(genId->at(iGen)) == abs(idGenB))
        evtWeight = 2;
    }
    
  }
  else  // If no gen information is available, only the BuToJPsiK channel can be useful
  {
    // Set id randomly for Bs - this should not be used
    if(selectionSubStrings[0].compare("BsToJPsiPhi") == 0)
    {
      idGenB = ((double)rand() / (RAND_MAX)) < 0.5 ? +531 : -531;
      std::cout << "EleMVASecondNtupleProducer::analyze():  W A R N I N G ! Trying to save tagging information for BsToJPsiPhi but\n";
      std::cout << "                                                        no generator information is available in the input files!\n";
      std::cout << "                                                        idGenB is filled randomly to " << idGenB << "!\n";
    }
    // Bu is a self-tagging channel - OK!
    if(selectionSubStrings[0].compare("BuToJPsiK") == 0)
    {
      for(auto iTrk: tracksFromB)
      {
        if(iTrk == tracksFromJPsi[0] || iTrk == tracksFromJPsi[1])
          continue;
        idGenB = trkCharge->at(iTrk) > 0 ? +521 : -521;
      }
      
    }
  }

  
  
  
  // Truth information for electron
  int iGenEle = -1;
  int idGenEle = 0;
  int genEleBMot = -1;
  
  // Variables needed to compute the tagging truth (stored into tagTruth)
  int chargeB = 0;
  int chargeCorr = 0;
  int chargeEle = eleCharge->at(iBestEle);
  
  // Use generator information if available
  if(has_gen)
  {
    // Truth information for electron
    iGenEle = GetClosestGen(elePt->at(iBestEle), eleEta->at(iBestEle), elePhi->at(iBestEle));
    if (iGenEle >= 0)
    {
      idGenEle = genId->at(iGenEle);
      genEleBMot = RecursiveLookForMother(iGenEle, allBHadrons);
      if(genEleBMot < 0)
        genEleBMot = -1;
    }
    if(selectionSubStrings[0].compare("BsToJPsiPhi") == 0)
    {
      if(abs(idGenEle) == 11)
        chargeCorr = GetGenLepBsChargeCorrelation(iGenEle, iGenB);
      if(!chargeCorr)
      {
        chargeCorr = GetBsChargeCorrelation(chargeEle, iGenB);
        chargeCorr*=2;
      }
    }
    else if(selectionSubStrings[0].compare("BuToJPsiK") == 0)
    {
      if(abs(idGenEle) == 11)
        chargeCorr = GetGenLepBuChargeCorrelation(iGenEle, iGenB);
      if(!chargeCorr)
      {
        chargeCorr = GetBuChargeCorrelation(chargeEle, iGenB);
        chargeCorr*=2;
      }
    }    
  }  
  else  // If no gen information is available, only the BuToJPsiK channel can be useful
  {
    // Bu is a self-tagging channel - OK!
    if(selectionSubStrings[0].compare("BuToJPsiK") == 0)
    {
      chargeB = 0;
      for(auto iTrack : tracksFromB)
      {
        chargeB += trkCharge->at(iTrack);
      }
      
      chargeCorr = chargeEle * chargeB;
    }    
  }
  
  int tagTruth;
  if(chargeCorr > 0)
    tagTruth = 1;
  else
    tagTruth = 0;
  
  std::cout << "Charge corr = " << chargeCorr << std::endl;
  std::cout << "Tag truth = " << tagTruth << std::endl;
  
  nAllTags++;
  if((chargeCorr > 0 && tagDecision < 0) || (chargeCorr < 0 && tagDecision > 0))
  {
    std::cout << "Different!\n\n";
    nMisTags++;
  }
  sumMisTagProb+=misTagProb;
  
  nGenB++; nGenB--;
  evtWeight++; evtWeight--;
  return true;
}


void OSElectronFlavorTagger::endJob() 
{
  // This runs after the event loop
  float misTagProbFromRatio = (float)nMisTags/(float)nAllTags;
  float misTagProbFromDNN = sumMisTagProb/(float)nAllTags;
  float dilutionFromRatio = 1. - 2. * misTagProbFromRatio;
  float dilutionFromDNN = 1. - 2. * misTagProbFromDNN; 
  float taggingEfficiency = (float)nAllTags/(float)nAllBEvts;
  float taggingPowerFromRatio = taggingEfficiency * dilutionFromRatio * dilutionFromRatio;
  float taggingPowerFromDNN = taggingEfficiency * dilutionFromDNN * dilutionFromDNN;
  std::cout << "tagging efficiency from ratio = " << taggingEfficiency << std::endl;
  std::cout << "misTagProb from dnn = " << misTagProbFromDNN << std::endl; 
  std::cout << "misTagProb from ratio = " << misTagProbFromRatio << std::endl; 
  std::cout << "dilution from dnn = " << dilutionFromDNN << std::endl; 
  std::cout << "dilution from ratio = " << dilutionFromRatio << std::endl;   
  std::cout << "tagging power from dnn = " << taggingPowerFromDNN << std::endl; 
  std::cout << "tagging power from ratio = " << taggingPowerFromRatio << std::endl; 
  return;
}
