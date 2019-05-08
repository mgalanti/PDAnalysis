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

  getUserParameter("tightSelection", tightSelection);

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
  
  // Set the generator-level B meson mass according to the particle looked for in the selection
  trueBMass = 0;
  if(evtSelection.substr(0,2).compare("Bs") == 0)
    trueBMass = constants::BsMass;
  else if(evtSelection.substr(0,2).compare("Bu") == 0)
    trueBMass = constants::BuMass;
  else
  {
    std::cout << "E R R O R ! Cannot determine the selected meson from the selection string!\n";
    std::cout << "            The selection string is: \"" << evtSelection << "\"\n";
    std::cout << "            Please fix the configuration file.\n";
    std::cout << "            Exiting...\n";
    exit(1);
  }
  std::cout << "EleMVASecondNtupleProducer::beginJob(): I N F O. trueBMass = " << trueBMass << std::endl;


  return;
}


void EleMVASecondNtupleProducer::book() 
{

  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility
  
  autoSavedObject = hEleBTrkDistance = Create2DHistogram<TH2D>("hEleBTrkDistance", "Distance in #Delta(R) and #Delta(p_{T})/p_{T} between electron and B tracks", 100, 0., 5., 100, 0., 2., "#Delta(R)", "2#times|p_{T,ele}-p_{T,Trk}|/(p_{T,ele}+p_{T,Trk})");

  return;
}


void EleMVASecondNtupleProducer::reset() {
// automatic reset
  autoReset();
  return;
}


bool EleMVASecondNtupleProducer::analyze(int entry, int event_file, int event_tot) 
{
  if (verbose) 
  {
    cout << " +++++++++++++++++++++++++++ " << endl;
    cout << "entry: "
         << entry << " " << event_file << " " << event_tot << endl;
    cout << "run: " <<   runNumber << " , "
         << "evt: " << eventNumber << endl;
  }
  else 
  {
//    if (!(event_file % 10000) || !(event_tot % 10000))
    if (!( event_tot % 10000 ) && event_tot)
      cout << event_file << " " << event_tot << endl;
  }
  if(entry == 0)
    checkBranches();

  // Main per-event analysis code goes here
  
  tWriter->Reset();
  
  // Check if the event passes the default selection defined in the cfg
  bool evtSelected = SelectEvent();
  
  if(!evtSelected)
    return false;
 
  // Check if the event passes the tight selection
  std::vector<std::pair<int, int> > selectedObjectsTight;
  bool tightEvent = SelectEvent(tightSelection, selectedObjectsTight);

  // If the event passes the tight selection, let's use those candidates instead of the loose ones
  if(tightEvent)
    selectedObjects = selectedObjectsTight;

  int iSelObject = 0;
  int iBestPV = -1;
  int iBestB = -1;
  int iBestEle = -1;
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
  
  std::vector <int> tracksFromB = tracksFromSV(iBestB);
  TLorentzVector pB = GetTLorentzVectorFromJPsiX(iBestB);
  
  // For monitoring purposes
  for(auto iTrk: tracksFromB)
  {
    float dR = deltaR(trkEta->at(iTrk), trkPhi->at(iTrk), eleGsfEta->at(iBestEle), eleGsfPhi->at(iBestEle));
    float dpTOverpT = 2*fabs(trkPt->at(iTrk) - eleGsfPt->at(iBestEle))/(trkPt->at(iTrk) + eleGsfPt->at(iBestEle));
    hEleBTrkDistance->Fill(dR, dpTOverpT);
  }
  
  // Signal-side variables
  (tWriter->tightEvent) = tightEvent?1:0;
  
  (tWriter->BPt) = pB.Pt();
  (tWriter->BEta) = pB.Eta();
  (tWriter->BPhi) = pB.Phi();
  (tWriter->BMass) = svtMass->at(iBestB);
  
  (tWriter->BLxy) = GetCt2D(pB, iBestB) / (trueBMass / pB.Pt());
  (tWriter->BCt2DBS) = GetCt2D(pB, iBestB, trueBMass);
  
  (tWriter->BCt2DPV) = GetCt2DPV(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt2DPVErr) = GetCt2DPVErr(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt2DPVSigmaUnit) = GetCt2DPV(pB, iBestB, iBestPV, trueBMass) / GetCt2DPVErr(pB, iBestB, iBestPV, trueBMass);

  (tWriter->BCt3DPV) = GetCt3DPV(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt3DPVErr) = GetCt3DPVErr(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt3DPVSigmaUnit) = GetCt3DPV(pB, iBestB, iBestPV, trueBMass) / GetCt3DPVErr(pB, iBestB, iBestPV, trueBMass);

  // Opposite-side variables
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
  
  autoSavedObject = cEleBTrkDistance = CreateCanvas("cEleBTrkDistance", "colz", false, false, false, hEleBTrkDistance);
  
  return;
}
