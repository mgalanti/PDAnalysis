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
  
  float minMass = 5.0;
  float maxMass = 5.5;
  int nMassBins = 250;
  
  autoSavedObject = hEleBTrkDistance = Create2DHistogram<TH2D>("hEleBTrkDistance", "Distance in #Delta(R) and #Delta(p_{T})/p_{T} between electron and B tracks", 100, 0., 5., 100, 0., 2., "#Delta(R)", "2#times|p_{T,ele}-p_{T,Trk}|/(p_{T,ele}+p_{T,Trk})");
  
  autoSavedObject = hWeightedBMass = Create1DHistogram<TH1D>("hWeightedBMass", "B mass (entries weighted by number of gen B in the event)", nMassBins, minMass, maxMass, "M(B) [GeV]", "Event weights");
  
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
  
  // HLT bit values
  bool jPsiMuHltBit = hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v) || hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v);
  bool jPsiTrkTrkHltBit = hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v);
  bool jPsiTrkHltBit = hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v);
  
  std::vector<int> tracksFromB = tracksFromSV(iBestB);
  TLorentzVector pB = GetTLorentzVectorFromJPsiX(iBestB);
  int iBestJPsi = (subVtxFromSV(iBestB)).at(0);
  std::vector<int> tracksFromJPsi = tracksFromSV(iBestJPsi);
  
  // For monitoring purposes
  for(auto iTrk: tracksFromB)
  {
    float dR = deltaR(trkEta->at(iTrk), trkPhi->at(iTrk), eleGsfEta->at(iBestEle), eleGsfPhi->at(iBestEle));
    float dpTOverpT = 2*fabs(trkPt->at(iTrk) - eleGsfPt->at(iBestEle))/(trkPt->at(iTrk) + eleGsfPt->at(iBestEle));
    hEleBTrkDistance->Fill(dR, dpTOverpT);
  }
  
  int iGenB = -1;
  int nGenB = -1;
  int idGenB = 0;
  int genBMixStatus = -1;
  int evtWeight = 1;
  // Use generator information if available
  if(has_gen)
  {
    std::vector<int> longLivedBHadrons = GetAllLongLivedBHadrons();
    iGenB = GetClosestGenInList(pB.Pt(), pB.Eta(), pB.Phi(), longLivedBHadrons, 0.4, 0.4);
    nGenB = longLivedBHadrons.size();
    // Do not consider events where the signal-side is not matched to a gen B
    if(iGenB < 0)
      return false;

    idGenB = genId->at(iGenB);

    // Do not consider events where the signal B is matched to a hadron different from the one looked for
    if(selectionSubStrings[0].compare("BsToJPsiPhi") == 0 && abs(idGenB) != 531)
    {
      std::cout << "Event rejected 1\n";
      return false;
    }
    if(selectionSubStrings[0].compare("BuToJPsiK") == 0 && abs(idGenB) != 521)
    {
      std::cout << "Event rejected 2\n";
      return false;
    }
    // Check whether the matched B has mixed or will mix
    genBMixStatus = GetMixStatus(iGenB);
    
    // FIXME: This is taken from Alberto's code, however I am not sure that it is the correct logic...
    if(genBMixStatus == 2)
    {
      std::cout << "Event rejected 3\n";
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
        evtWeight++; // FIXME: was evtWeight = 2 in original code. Check which is correct.
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
  
  // General event variables
  (tWriter->evtNumber) = event_tot;
  (tWriter->evtWeight) = evtWeight;
  (tWriter->tightEvent) = tightEvent;

  (tWriter->nGenB) = nGenB;
  
  (tWriter->JPsiMuHltBit) = jPsiMuHltBit;
  (tWriter->JPsiTrkTrkHltBit) = jPsiTrkTrkHltBit;
  (tWriter->JPsiTrkHltBit) = jPsiTrkHltBit;  
  
  (tWriter->iPV) = iBestPV;
    
  // Signal-side variables
  (tWriter->BPt) = pB.Pt();
  (tWriter->BEta) = pB.Eta();
  (tWriter->BPhi) = pB.Phi();
  (tWriter->BMass) = svtMass->at(iBestB);
  
  hWeightedBMass->Fill(svtMass->at(iBestB), evtWeight);
  
  (tWriter->BLxy) = GetCt2D(pB, iBestB) / (trueBMass / pB.Pt());
  (tWriter->BCt2DBS) = GetCt2D(pB, iBestB, trueBMass);
  
  (tWriter->BCt2DPV) = GetCt2DPV(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt2DPVErr) = GetCt2DPVErr(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt2DPVSigmaUnit) = GetCt2DPV(pB, iBestB, iBestPV, trueBMass) / GetCt2DPVErr(pB, iBestB, iBestPV, trueBMass);

  (tWriter->BCt3DPV) = GetCt3DPV(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt3DPVErr) = GetCt3DPVErr(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt3DPVSigmaUnit) = GetCt3DPV(pB, iBestB, iBestPV, trueBMass) / GetCt3DPVErr(pB, iBestB, iBestPV, trueBMass);

  (tWriter->BiSV) = iBestB;
  
  (tWriter->BidGen) = idGenB;
  
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

  autoSavedObject = cWeightedBMass = CreateCanvas("cWeightedBMass", 0, 21, 1, false, false, hWeightedBMass);
  
  return;
}
